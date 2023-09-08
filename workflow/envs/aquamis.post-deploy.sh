#!/usr/bin/env bash
set -Eeu

# Repo URL
aquamis_repo="https://gitlab.com/bfr_bioinformatics/AQUAMIS.git"

# Commit hash to use
commit="0290eb2acafd697f631b8b3a6135d2710519b575"

# Local directory to save the Repo
local_dir="$HOME/.nrw-geuebt/aquamis"

# if already exists, wipe it clean and redo clone
[ -d "$local_dir" ] && rm -rf "$local_dir"

# clone and checkout repo
echo "Cloning AQUAMIS and checking out stable commit"
git clone -q "$aquamis_repo" "$local_dir"
cd "$local_dir"
git checkout "$commit"

# Patch wrapper and Snakefile
# echo "Applying patches"
# patch -s --directory="$local_dir" --strip=1 << END
# diff --unified --recursive --no-dereference chewieSnake-orig/chewieSnake_join.py chewieSnake/chewieSnake_join.py
# --- chewieSnake-orig/chewieSnake_join.py	2023-05-30 10:12:16.855136923 +0200
# +++ chewieSnake/chewieSnake_join.py	2023-05-30 10:14:42.000000000 +0200
# @@ -147,7 +147,7 @@
     # parser.add_argument('--clustering_method', help = 'The agglomeration method to be used for hierarchical clustering. This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC); default: single',
                         # default = "single", required = False)
     # parser.add_argument('--distance_threshold', help = 'A single distance threshold for the exctraction of sub-trees; default: 10',
# -                        default = 10, required = False)
# +                        default = 10, type=int, required = False)
     # parser.add_argument('--cluster_names', help = 'A file with potential names for cluster names, one name per line, default: repo_path/scripts/cluster_names_reservoir.txt',
                         # default = os.path.join(repo_path, "scripts", "cluster_names_reservoir.txt"), required = False)
     # parser.add_argument('--subcluster_thresholds', help = 'A list of distance thresholds for subclustering; default: [3]',
# diff --unified --recursive --no-dereference chewieSnake-orig/chewieSnake_join.smk chewieSnake/chewieSnake_join.smk
# --- chewieSnake-orig/chewieSnake_join.smk	2023-05-30 10:12:16.451155752 +0200
# +++ chewieSnake/chewieSnake_join.smk	2023-05-30 10:14:52.000000000 +0200
# @@ -76,7 +76,7 @@
         # "merged_db/sample_cluster_information.tsv", # from rule clustering
         # "merged_db/samples2clusters.csv" # result from checkpoint
         # "merged_db/listofreports.txt", # result from checkpoint
# -        "merged_db/report.html"
# +        # "merged_db/report.html"
 
 # rule all_reportonly:
     # input:

# END

# update env with aquamis.yaml
if [[ $(which mamba) == "" ]]; then
  conda env update --file "$local_dir/envs/auamis.yaml" -p $CONDA_PREFIX &> "$local_dir/deploy_log.txt"
  # conda update -y -p $CONDA_PREFIX "python>=3.9" &>> "$local_dir/deploy_log.txt"
else
  mamba env update --file "$local_dir/envs/aquamis.yaml" -p $CONDA_PREFIX &> "$local_dir/deploy_log.txt"
  # mamba update -y -p $CONDA_PREFIX "python>=3.9" &>> "$local_dir/deploy_log.txt"
fi

# Install databases - using aquamis script
bash "$local_dir/scripts/aquamis_setup.sh" --databases &>> "$local_dir/deploy_log.txt"

# Finalize BUSCO install - modified from scirpt to work with hashed env name
download_file() {
  # Init
  download_success=false
  download_status=false
  download_hash=''
  
  # Download
  if [[ -f $local_archive ]] && [[ $arg_force == false ]]; then
    echo "The file $local_archive already exists. Skipping download."
  else
    echo "Downloading $remote_archive to $local_archive"
    wget --ca-certificate=$local_dir/resources/seafile-bfr-berlin_certificate.pem --output-document $local_archive $remote_archive && download_status=true
    [[ "$download_status" = true ]] && [[ -s $local_archive ]] && download_hash=$(openssl dgst -r -sha256 $local_archive) && download_success=true
  fi
  
  # Verify Integrity of Download
  [[ -n $download_hash ]] && echo "$download_hash" >> $hashfile

  # Unpack Downloads
  [[ "$download_success" == true ]] && echo "Extracting archive $(basename $local_archive) to $extraction_directory"
  [[ "$download_success" == true ]] && tar -xzv -f $local_archive -C $extraction_directory
}

mkdir -p "$local_dir/download"

remote_archive="https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/augustus.tar.gz" # 139MB or https://databay.bfrlab.de/f/fa13c9eb2625477eb729/?dl=1
local_archive="$local_dir/download/augustus.tar.gz"
hashfile="$local_dir/download/reference_db.sha256"
extraction_directory="$CONDA_PREFIX/lib/python3.7/site-packages/quast_libs/"
download_file  &>> "$local_dir/deploy_log.txt"

remote_archive="https://gitlab.bfr.berlin/bfr_bioinformatics/aquamis_databases/-/raw/main/bacteria.tar.gz" # 8.8MB or https://busco.ezlab.org/v2/datasets/bacteria_odb9.tar.gz
local_archive="$local_dir/download/bacteria.tar.gz"
hashfile="$local_dir/download/reference_db.sha256"
extraction_directory="$CONDA_PREFIX/lib/python3.7/site-packages/quast_libs/busco/"
download_file  &>> "$local_dir/deploy_log.txt"

cat "$local_dir/download/reference_db.sha256" >> "$local_dir/reference_db/reference_db.sha256"
rm -r "$local_dir/download"
