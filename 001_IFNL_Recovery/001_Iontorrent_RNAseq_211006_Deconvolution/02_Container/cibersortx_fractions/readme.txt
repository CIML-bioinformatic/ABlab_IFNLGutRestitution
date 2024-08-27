
This image contains:

 - Cibersortx Fractions


# ######################
     COMPILE THE IMAGE
# ######################

docker build -t ablab_iflnrecovery_cibersort /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/spinelli/ciml-ablab/Project/001_IFNL_Recovery/001_Iontorrent_RNAseq_211006_Deconvolution/02_Container/cibersortx_fractions/

# #############################################
# Start the container
# #############################################

docker run -v {dir_path}:/src/data -v {dir_path}:/src/outdir cibersortx/fractions [Options]



