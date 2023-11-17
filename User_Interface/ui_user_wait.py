""" A Small Continue/Halt dialog box that utilizes Return Objects for the result
@author: Ben C
"""
import tkinter as tk
import mcn_status as mcs


class Waiter:
    """ A basic pause dialogue box """
    def __init__(self, root: tk.Tk, info):
        """ Creates a small prompt that pauses execution until dismissed

        :param root: a tkinter root
        :param info: the prompt
        """
        self.core = root
        frame = tk.Frame(self.core)
        frame.winfo_toplevel().title("User Feedback Required")
        frame.pack(side=tk.LEFT, fill=tk.BOTH, padx=10)
        self.ret_val = mcs.RetObj.incomplete("_U", mcs.V_FATAL, "Unhandled GUI Exit Value")
        prompt = tk.Label(frame, text=info)
        prompt.grid(row=0, column=0, columnspan=2)
        continue_button = tk.Button(frame, text="Continue", command=lambda: self.cont())
        continue_button.grid(row=1, column=0)
        halt_button = tk.Button(frame, text="Halt", command=lambda: self.halt())
        halt_button.grid(row=1, column=1)

    def cont(self):
        """ Sets the return to a complete Return Object, then closes

        :return: None
        """
        self.ret_val = mcs.RetObj.complete("User told MCN to continue")
        self.core.destroy()

    def halt(self):
        """ Sets the return to an incomplete Return Object, then closes

        :return: None
        """
        self.ret_val = mcs.RetObj.incomplete("_U", mcs.V_FATAL, "User issued Halt")
        self.core.destroy()

    def run(self) -> mcs.RetObj:
        """ Brings up the poppup

        :return: A return object (complete if "continue", incomplete otherwise)
        """
        self.core.mainloop()
        return self.ret_val


if __name__ == '__main__':
    _root = tk.Tk()
    _info = "Custom Message to Prompt User"
    w = Waiter(_root, _info)
    w.run()
