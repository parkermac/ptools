"""
Code to tesk tkinter.

Must be run BEFORE invoging matplotlib, so on my mac you can
run it in ipython, but not ipy='ipython --pylab'

From:
http://interactivepython.org/courselib/static/thinkcspy/
GUIandEventDrivenProgramming/01_intro_gui.html#tkinter

"""

import tkinter as tk
from tkinter import ttk

# Create the application window
window = tk.Tk()

# Create the user interface
my_label = ttk.Label(window, text="Hello World!")
my_label.grid(row=1, column=1)

# Start the GUI event loop
window.mainloop()