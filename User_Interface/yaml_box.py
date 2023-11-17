import yaml
from gui_constants import *
import tkinter as tk
from io import StringIO
from typing import Optional


class YamlBox:
    def __init__(self, core, w=18, h=5, stretchy=False):
        self.core = core
        self.frame = tk.Frame(self.core)
        # self.frame.pack(**TB33)
        self.container = tk.Frame(self.frame)
        self.container.pack(**TB33, expand=stretchy)
        self.scrolly = tk.Scrollbar(self.container)
        self.scrollx = tk.Scrollbar(self.container, orient='horizontal')
        self.scrolly.pack(**PACK_SCROLL)
        self.scrollx.pack(side=tk.BOTTOM, fill=tk.X)
        self.textbox = tk.Text(self.container,
                               yscrollcommand=self.scrolly.set,
                               xscrollcommand=self.scrollx.set,
                               width=w,
                               height=h,
                               wrap=tk.NONE)
        self.textbox.pack(**LB33, expand=stretchy)
        self.scrolly.config(command=self.textbox.yview)
        self.scrollx.config(command=self.textbox.xview)
        self.notif_frame = tk.Frame(self.frame)
        self.notif_frame.pack(**TB33, expand=stretchy)
        self.notification = tk.Label(self.notif_frame, text='--')
        self.notification.pack(**BOTH33, expand=stretchy)

    def get(self, start, end) -> Optional[dict]:
        text = self.textbox.get(start, end)
        try:
            doc = yaml.safe_load(text)
        except yaml.YAMLError as ye:
            mark = ye.problem_mark  # noqa # it's fine, it's there
            line = mark.line + 1
            column = mark.column + 1
            self.notification['text'] = f'YAML Error: L{line} C{column}'
            return None
        else:
            self.notification['text'] = '--'
            return doc

    def delete(self, start, end):
        return self.textbox.delete(start, end)

    def insert(self, where, what):
        stream = StringIO()
        yaml.safe_dump(what, stream)
        self.textbox.insert(where, stream.getvalue())
        del stream

    def pack(self, *args, **kwargs):
        return self.frame.pack(*args, **kwargs)

    def grid(self, *args, **kwargs):
        return self.frame.grid(*args, **kwargs)
