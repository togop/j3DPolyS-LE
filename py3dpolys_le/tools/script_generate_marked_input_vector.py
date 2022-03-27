import boundary_files as bdf
import drawSvg as draw


def draw_input_vector_svg(vector_of_genomic_locati, filepath_chromosome_sizes, output_svg_filename,
                          chromosome_of_interest="chrX",
                          resolution=10000,
                          height_width_ratio=24):

    """
    :param vector_of_genomic_locati: locations that will show vertical lines on the svg
    :param filepath_chromosome_sizes: chromosomesizes.txt file used to get the chromosome's length
    :param chromosome_of_interest: which chromosome to pick from the file
    :param resolution: genomic position will be divided by this factor
    :param height_width_ratio: the aspect ratio of the resulting svg, 24 times wider than tall looks good
    :return: doesn't return anything, but writes a .svg in the working directory
    """

    # get data
    vector_of_genomic_locati = [round(location / resolution) for location in vector_of_genomic_locati]
    print(vector_of_genomic_locati)

    chr_sizes = bdf.read_chromosome_sizes(filepath_chromosome_sizes)
    chr_size = round(int(chr_sizes[chromosome_of_interest]) / resolution)
    print(chr_size)

    width = chr_size
    height = round(chr_size / height_width_ratio)

    svg_drawing = draw.Drawing(width, height, displayInline=False)

    # draw chromosome from start to finish
    svg_drawing.append(draw.Lines(0, (height / 2), width, (height / 2), close=False, stroke='black'))

    # draw perpendicular line for each boundary
    for boundary in vector_of_genomic_locati:
        svg_drawing.append(draw.Lines(boundary, 0, boundary, height, close=False, stroke='black'))

    svg_drawing.saveSvg(output_svg_filename + ".svg")


def generate_vectorSVG_from_boundaryfile(filepath_boundaryfile, filepath_chromosome_sizes,
                                         chromosome_of_interest="chrX",
                                         resolution=10000,
                                         height_width_ratio=24):

    boundary_arr = bdf.read(filepath_boundaryfile)
    boundaries = [int(boundary[1]) for boundary in boundary_arr]

    filename = str(filepath_boundaryfile.split(sep="/")[-1].split(".")[0])

    draw_input_vector_svg(boundaries, filepath_chromosome_sizes, chromosome_of_interest, filename, resolution, height_width_ratio)


if __name__ == "__main__":

    # parameters
    filepath_boundaryfile = "../../py3dpolys_le/data/ce/boundary_sites/boundaries_chrX_manual_detection.csv"
    filepath_outfile = "../../py3dpolys_le/data/simulations/chip_seqs/boundary_mnl_bid_asym.chip.out"

    filepath_chromosome_sizes = "../../py3dpolys_le/data/ce/chromosome_sizes_chrIVXnotation.tsv"

    generate_vectorSVG_from_boundaryfile(filepath_boundaryfile, filepath_chromosome_sizes)




