import numpy as np
import pandas as pd

class SourceSinkIn():
    """ class for *.prop or other similar formats"""
    def __init__(self, filename=None, number_of_groups=2, ele_groups=[]):
        self.n_group = number_of_groups  # 0: source; 1: sink
        if filename is not None:
            """Read a usgs data file and initialize from its content """
            self.source_file = filename
            self.np_group = []
            self.ip_group = []

            with open(self.source_file) as fin:
                for k in range(0, self.n_group):
                    self.np_group.append(int(fin.readline().split()[0]))
                    print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                    self.ip_group.append(np.empty((self.np_group[k]), dtype=int))
                    for i in range(0, self.np_group[k]):
                        self.ip_group[k][i] = int(fin.readline())
                    fin.readline()  # empty line between groups
                    if self.np_group[k] > 0:
                        print("p first: " + str(self.ip_group[k][0]))
                        print("p last: " + str(self.ip_group[k][-1]))
        else:
            self.np_group = [len(x) for x in ele_groups]
            self.ip_group = [np.array(x) for x in ele_groups]

    def writer(self, filename=None):
        if filename is None:
            filename = self.source_file

        with open(filename, 'w') as fout:
            for k in range(0, self.n_group):
                print("Points in Group " + str(k+1) + ": " + str(self.np_group[k]))
                fout.write(f"{self.np_group[k]}\n")
                for i in range(0, self.np_group[k]):
                    fout.write(f"{self.ip_group[k][i]}\n")
                fout.write("\n")  # empty line

def relocate_sources(old_source_sink_in:SourceSinkIn, old_vsource:np.ndarray, outdir=None, relocate_map = None):

    # make a dataframe from vsource to facilitate column subsetting and reordering
    df = pd.DataFrame(data=old_vsource, columns=['time'] + [str(x) for x in old_source_sink_in.ip_group[0]])

    # read relocate_map
    if relocate_map is not None:
        eleids = relocate_map[:, 0]; new2old_sources = relocate_map[:, 1]
    else:
        raise ValueError('relocate_map is not provided')

    # Assemble new source_sink.in
    source_sink_in = SourceSinkIn(filename=None, number_of_groups=2, ele_groups=[eleids.tolist(),[]])
    source_sink_in.writer(f'{outdir}/source_sink.in')

    # Rearrange vsource.th and msource.th
    map_dict = {str(k): str(v) for k, v in np.c_[old_source_sink_in.ip_group[0][new2old_sources], eleids]}
    map_dict = {**{'time': 'time'}, **map_dict}

    # Check if all columns in mapping exist in the dataframe
    assert set(map_dict.keys()).issubset(set(df.columns)), "Some columns in the mapping don't exist in the DataFrame"

    # Subset and rename the columns based on the mapping
    df_subset = df[list(map_dict.keys())].rename(columns=map_dict)

    # Save the subset dataframe to a new CSV file
    df_subset.to_csv(f'{outdir}/vsource.th', index=False, header=False, sep=' ')

    # msource
    msource = np.c_[
        np.r_[df['time'].iloc[0], df['time'].iloc[-1]],
        np.ones((2, len(eleids)), dtype=int) * -9999,
        np.zeros((2, len(eleids)), dtype=int)
    ]
    np.savetxt(f'{outdir}/msource.th', msource, fmt='%d', delimiter=' ')

    pass

if __name__ == "__main__":
    old_source_sink_in = SourceSinkIn(filename='../original_source_sink/source_sink.in')
    old_vsource = np.loadtxt('../original_source_sink/vsource.th')
    relocate_map = np.loadtxt(f'./relocate_map.txt', dtype=int)
    outdir = './'

    relocate_sources(
        old_source_sink_in=old_source_sink_in,
        old_vsource=old_vsource,
        outdir=outdir,
        relocate_map=relocate_map
    )
