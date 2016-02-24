from imports import *

dirname = "/Users/omidi/Dropbox/AFM Data/03-09-2014_HeightData"
filenames = listdir(dirname)


def load_file(filename):
    position = []
    height = []
    with open(path.join(dirname, filename), "rU") as inf:
        for record in DictReader(inf, dialect='excel', delimiter=","):
            position.append(float(record["X"]))
            height.append(float(record["Y"]))
    return np.array(position), np.array(height)


def return_all_files_for_cell(cellname):
    return [f for f in filenames if re.search("_%s.csv$" % cellname, f)]
