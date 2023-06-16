library("rhdf5")
library("optparse")
library("combinat")
# library("imager")
library(stringr)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input file."),
  make_option(c("-o", "--output"), type = "character", help = "Output file's prefix.")
)

parseobj <- OptionParser(option_list=option_list)
opt <- parse_args(parseobj)
print(opt)

# just for testing purpose
###### Using already merged ##########
input_file <- if (!is.null(opt$input)) opt$input else "/Users/todor/data/3dpolys_le/simulations/demo/out-Nlef225-km0.0027-bd0_im-z/r2.14/hic3d_003_10k_cool.hdf5"
output_prefix <- if (!is.null(opt$output)) opt$output else "/Users/todor/data/3dpolys_le/simulations/demo/out-Nlef225-km0.0027-bd0_im-z/r2.14/hic3d/hic3d_003_10k"

if (!dir.exists(dirname(output_prefix))){
  dir.create(dirname(output_prefix))
}

# hic3d_fid <- H5Fopen(input_file)

## hic3d_h5 <- H5Dopen(hic3d_fid, "hic3d_cool")
hic3_cool <- h5read(file = input_file, name = "hic3d_cool")
# hic3d_bins <- h5read(file = input_file, name = "bins")

# hic3_cool <- hic3_cool
# hic3d_bins <- t(hic3d_bins)

#bins_df <- data.frame(hic3d_bins)  # s, x, y, c)
#colnames(bins_df) <- c("chrom", "start", "end")
#bins_filename <- paste0(output_prefix, '_bins.tsv')
#write.table(bins_df, file=bins_filename, sep='\t', row.names = FALSE, quote=FALSE)

hic3_N <- max(hic3_cool[, 3]) # , max(hic3_cool[,2]), max(hic3_cool[,3]))  # dim(hic3d_bins)[1]
max_c <- max(hic3_cool[, 4])
col_pal <- rev(grey(seq(0, 1, length = max_c)))
zplots_df <- data.frame(matrix(ncol = 4, nrow = 0))  # s, x, y, c)
colnames(zplots_df) <- c("s", "x", "y", "c")
zplot_hic <- matrix(0, nrow=hic3_N, ncol=hic3_N)
# fileter: 0,0,0,0
hic3_cool <- hic3_cool[hic3_cool[, 4] > 0, ]
for (si in (1:hic3_N)){
  # si_zplot_h5 <- H5Dopen(paste0(output_prefix, '_s', si, '.cool'), "hic3d_cool")
  zplot <- matrix(0, nrow=hic3_N, ncol=hic3_N)
  # image(zplot, useRaster=FALSE, axes=FALSE, col = grey(seq(0, 1, length = 256)))
  for (xyz in permn(3)){
    xi <- xyz[1]
    yi <- xyz[2]
    zi <- xyz[3]
    s_ids <- hic3_cool[, zi] == si
    if (any(s_ids)){
      x <- hic3_cool[s_ids, xi]
      y <- hic3_cool[s_ids, yi]
      c <- hic3_cool[s_ids, 4]
      for (i in seq(length(x))){
        if ((x[i] <= hic3_N) && (y[i] <= hic3_N)) {
          zplot[x[i], y[i]] <- c[i]
        } else {
          cat('zplot out of bound x,y: ', x[i], y[i], '\n')
        }
      }
    }

    xyids <- which(zplot!=0,arr.ind = T)
    s <- rep(si, length(xyids))
    x <- xyids[,1]
    y <- xyids[,2]
    c <- zplot[xyids]

    #sub_zplot_df <- data.frame(s, x, y, c)
    #zplots_df <- rbind(zplots_df, sub_zplot_df)
  }
  filename <- paste0(output_prefix, '_s', str_pad(si, 3, pad = "0"), '.tif')
  #if (!isSymmetric(zplot)){
  #  cat('Not simmetric ', filename, '\n')
  #}
  cat('save image ', filename, '\n')
#  imager::save.image(zplot, filename)
  tiff(filename)
  image(zplot, useRaster=FALSE, axes=FALSE, col = col_pal)
  dev.off()

  zplot_hic <- zplot_hic + zplot
}

tiff(paste0(output_prefix, '_zhic.tif'))
image(zplot_hic, useRaster=FALSE, axes=FALSE, col = col_pal)
dev.off()

write.table(zplot_hic, file=paste0(output_prefix, '_zhic.tsv'), sep = "\t", row.names=FALSE, col.names=FALSE)

H5close()

# zplots_sorted_df <- zplots_df[order(zplots_df$s, zplots_df$x, zplots_df$y),]
#
# filename <- paste0(output_prefix, '_contacts.tsv')
# cat('Saving as .tsv ...', filename, '\n')
# write.table(zplots_sorted_df, file=filename, sep='\t', row.names = FALSE, quote=FALSE)
#
# cat('Saving as .tsv.gz ...\n')
# gz1 <- gzfile(paste0(filename,".gz"), "w")
# write.csv(zplots_sorted_df, gz1)
# close(gz1)
# cat('Sinished ',filename,".gz")

