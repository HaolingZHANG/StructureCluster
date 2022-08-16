from datetime import datetime


class Monitor(object):

    def __init__(self):
        """
        Initialize the monitor to identify the task progress.

        Example
            >>> from dsw import Monitor
            >>> monitor = Monitor()
            >>> monitor.output(current_state=1, total_state=10)
            \r|███                 | 10% ( 1/10) wait 0000:00:00.
            >>> monitor.output(current_state=5, total_state=10)
            \r|███████████         | 50% ( 5/10) wait 0000:00:00.
            >>> monitor.output(current_state=10, total_state=10)
            \r|████████████████████|100% (10/10) used 0000:00:00.
        """
        self.last_time = None

    def output(self, current_state, total_state, extra=None):
        """
        Output the current state of process.

        :param current_state: current state of process.
        :type current_state: int

        :param total_state: total state of process.
        :type total_state: int

        :param extra: extra vision information if required.
        :type extra: dict
        """
        if self.last_time is None:
            self.last_time = datetime.now()

        if current_state == 0:
            return

        position = int(current_state / total_state * 100)

        string = "|"

        for index in range(0, 100, 5):
            if position >= index:
                string += "█"
            else:
                string += " "

        string += "|"

        pass_time = (datetime.now() - self.last_time).total_seconds()
        wait_time = int(pass_time * (total_state - current_state) / current_state)

        string += " " * (3 - len(str(position))) + str(position) + "% ("

        string += " " * (len(str(total_state)) - len(str(current_state))) + str(current_state) + "/" + str(total_state)

        if current_state < total_state:
            minute, second = divmod(wait_time, 60)
            hour, minute = divmod(minute, 60)
            string += ") wait " + "%04d:%02d:%02d" % (hour, minute, second)
        else:
            minute, second = divmod(pass_time, 60)
            hour, minute = divmod(minute, 60)
            string += ") used " + "%04d:%02d:%02d" % (hour, minute, second)

        if extra is not None:
            string += " " + str(extra).replace("\'", "").replace("{", "(").replace("}", ")") + "."
        else:
            string += "."

        print("\r" + string, end="", flush=True)

        if current_state >= total_state:
            self.last_time = None
            print()
