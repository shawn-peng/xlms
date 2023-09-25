import matplotlib.pyplot as plt
import multiprocessing as mp
import matplotlib

# matplotlib.use('TkAgg', force=True)
matplotlib.use('Qt5Agg')
print("Switched to:", matplotlib.get_backend())



class ProcessPlotter(object):
    def __init__(self, plot_func):
        self.x = []
        self.y = []
        self.plot_func = plot_func
        # self.fig.show
        self.static_pos_args = {}
        self.static_kw_args = {}

    def terminate(self):
        plt.close('all')

    def call_back(self):
        count = 0
        while self.pipe.poll():
            command = self.pipe.recv()
            print('command', command)
            if command is None:
                self.terminate()
                return False
            elif command[0] == 'set_static_pos_arg':
                pos, arg = command[1:]
                self.static_pos_args[pos] = arg
            count += 1
            # self.x.append(command[0])
            # self.y.append(command[1])
            # self.ax.plot(self.x, self.y, 'ro')
        # print(count, 'commands read')
        if not count:
            return True
        cmd, *args = command
        args = list(args)
        # args = list(command[1:])
        if cmd == 'plot':
            for i, arg in self.static_pos_args.items():
                if args[i] is None:
                    args[i] = arg
            self.plot_func(*args)
            self.fig.canvas.draw()
            self.fig.canvas.flush_events()
        return True

    def __call__(self, pipe):
        print('starting plotter...')

        print('pipe', pipe)
        self.pipe = pipe
        self.fig = plt.figure(figsize=[18, 9])
        # self.fig, self.ax = plt.subplots()
        plt.ioff()
        timer = self.fig.canvas.new_timer(interval=10000)
        timer.add_callback(self.call_back)
        timer.start()

        print('...done')
        plt.show()
        print('plt.show finished')


class NBPlot(object):
    def __init__(self):
        self.plot_pipe, plotter_pipe = mp.Pipe()
        self.plotter = ProcessPlotter()
        self.plot_process = mp.Process(
            target=self.plotter, args=(plotter_pipe,), daemon=True)
        self.plot_process.start()

    def plot(self, finished=False):
        send = self.plot_pipe.send
        if finished:
            send(None)
        else:
            data = 1
            send(data)
