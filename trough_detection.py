from imports import *
from toolbox import *


def plot_one_cell_in_time(cell_name):
    data = map(load_file, return_all_files_for_cell(cell_name))
    orders_by_length = np.argsort(([max(d[0]) for d in data]))[::-1]
    longest_cell = orders_by_length[0]
    plt.plot(data[longest_cell][0], data[longest_cell][1] / max(data[longest_cell][1]))
    total_pad = 0
    for t in xrange(1, len(data)):
        t_prev = orders_by_length[t-1]
        t_curr = orders_by_length[t]
        x = data[t_curr][1]
        y = data[t_prev][1]
        conv = np.convolve(y, x[::-1])
        left_pad = np.argmax(conv) - len(x) + 1
        if left_pad<0:
            left_pad = 0
        total_pad += data[t_prev][0][left_pad]
        plt.plot(data[t_curr][0] + total_pad, (x / max(x)))
    plt.show()


plot_one_cell_in_time('1.1.1')