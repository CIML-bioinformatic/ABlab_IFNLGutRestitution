
This image contains:

 - Rstudio
 - R 4.1.1

# ######################
     COMPILE THE IMAGE
# ######################

docker build -t ablab_iflnrecovery_rstudio /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-ablab/Project/001_IFNL_Recovery/001_Iontorrent_RNAseq_211006_Deconvolution/02_Container/rstudio/

# #############################################
# Start the container
# #############################################

docker run --name ablab_iflnrecovery_rstudio -d -p 9292:8787 -v /mnt:/mnt -e PASSWORD= -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g)  ablab_iflnrecovery_rstudio



