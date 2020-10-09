# For each graded_cleaned_scrubbed_clipped_whatever...
# run variability_selection.spreadsheet_maker

from wuvars.analysis.variability_selection import spreadsheet_maker

def write_summary_spreadsheet(filename, output):
    # takes in a dataset from a given filename
    # makes an intermediate summary spreadsheet
    # writes that to file
    # returns nothing

    pass


if __name__ == "__main__":
    # do some globbing? didn't we do something like this before?
    # grab all of the things and then 
    input_filenames = []
    output_filenames = []

    for in_file, out_file in zip(input_filenames, output_filenames):

        write_summary_spreadsheet(in_file, out_file)