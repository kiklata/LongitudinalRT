## use cellbender env
# python 3.7
# pip install --no-cache-dir -U git+https://github.com/broadinstitute/CellBender.git@7fd0dac8fe5c37e705cdd50fa5767064f8f4b980
## google colab failed due to limited memory (13GB) of colab, run it in local. see https://github.com/broadinstitute/CellBender/issues/266

## just download all file required in local PC, run in cellbender env, succeed                                               
conda run -n cellbender cellbender remove-background --input raw_feature_bc_matrix.h5 \
                                                     --output cellbender_feature_bc_matrix.h5 \
                                                     --cuda

## P1013S2 UMAP display a poor perform, might not use it