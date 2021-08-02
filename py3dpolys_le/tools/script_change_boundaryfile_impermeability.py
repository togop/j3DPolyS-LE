import boundary_files as bdf


if __name__ == "__main__":

    filepath_boundaryfile_input = "/tmp/boundaries_chrX_mnl.csv"
    filepath_boundaryfile_output = "/tmp/boundaries_chrX_mnl.imp_0.95.csv"

    array_bdf = bdf.read(filepath_boundaryfile_input)
    array_bdf_mod = bdf.blanket_change_impermeability(array_bdf, 0.95)
    bdf.write(array_bdf_mod, filepath_boundaryfile_output)
