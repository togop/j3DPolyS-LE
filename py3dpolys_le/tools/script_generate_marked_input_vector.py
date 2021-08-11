import boundary_files as bdf
import drawSvg as draw

if __name__ == "__main__":

    # parameters
    filepath_boundaryfile = "../../py3dpolys_le/data/ce/boundary_sites/boundaries_chrX_manual_detection.csv"
    filepath_chromosome_sizes = "../../py3dpolys_le/data/ce/chromosome_sizes_chrIVXnotation.tsv"

    chromosome_of_interest = "chrX"
    resolution = 10000  # genomic position will be divided by this factor
    height_width_ratio = 24  # the aspect ratio of the resulting svg, 1:24 looks good

    # get data
    filename = str(filepath_boundaryfile.split(sep="/")[-1].split(".")[0])

    boundary_arr = bdf.read(filepath_boundaryfile)
    boundaries = [int(boundary[1]) for boundary in boundary_arr]
    boundaries = [round(boundary / resolution) for boundary in boundaries]
    print(boundaries)

    chr_sizes = bdf.read_chromosome_sizes(filepath_chromosome_sizes)
    chr_size = round(int(chr_sizes[chromosome_of_interest]) / resolution)
    print(chr_size)

    #
    width = chr_size
    height = round(chr_size / height_width_ratio)

    svg_drawing = draw.Drawing(width, height, displayInline=False)

    # draw chromosome from start to finish
    svg_drawing.append(draw.Lines(0, (height / 2), width, (height / 2), close=False, stroke='black'))

    # draw perpendicular line for each boundary
    for boundary in boundaries:
        svg_drawing.append(draw.Lines(boundary, 0, boundary, height, close=False, stroke='black'))

    svg_drawing.saveSvg(filename + ".svg")
