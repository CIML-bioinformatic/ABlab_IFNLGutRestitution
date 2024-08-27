
This image contains:

 - Rstudio
 - R 4.1.1

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t ablab_iflnrecovery_cduc /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-ablab/Project/001_IFNL_Recovery/002_IBD_CD_UC/02_Container/rstudio/

# #############################################
# Start the container
# #############################################

docker run --name ablab_iflnrecovery_cduc -d -p 9393:8787 -v /mnt:/mnt -e PASSWORD= -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  ablab_iflnrecovery_cduc



