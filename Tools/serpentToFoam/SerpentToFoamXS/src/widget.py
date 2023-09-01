"""
SerpentToFoamXS

Generic widgets definitions

Author: Thomas Guilbaud, EPFL/Transmutex SA
Last Update: 18/03/2022
"""

# Imports
import tkinter as tk
from tkinter import filedialog
import tkinter.ttk as ttk
from PIL import Image, ImageTk

from . import utils
from . import image

# Image
import base64
import io


#------------------------------------------------------------------------------*
# Filename input
# File Reader
def setFilenameSelected(filename, event=None) -> None:
    filename.set(tk.filedialog.askopenfilename())

def fileSelection(frame, filename, labelText: str, row: int, column: int=2) -> int:
    """
    Widget to select a file using a popup window.
    """
    # Avoid empty filename
    if (filename.get() == ""):
        filename.set("Open")

    ttk.Label(frame, text=labelText).grid(row=row, column=0, sticky="W")
    style = ttk.Style()
    style.configure('File.TButton', anchor=tk.E)
    ttk.Button(
        frame,
        textvariable=filename,
        width=25,
        style="File.TButton",
        command=lambda: setFilenameSelected(filename)
    ).grid(row=row, column=column, sticky="W")

    return(row+1)


#------------------------------------------------------------------------------*
# Widget Line

def line(masterFrame, labelText: str, variable, default, row: int, column: int=2, unit: str=""):
    """
    Standard input line. Choice between regular text entry or boolean checkbox
    using respectively tk.StringVar and tk.IntVar.
    """
    label = ttk.Label(masterFrame, text=labelText)
    label.grid(row=row, column=0, sticky="W")
    if (utils.isIntVar(variable)):
        entry = ttk.Checkbutton(masterFrame, variable=variable)
    elif (utils.isStringVar(variable)):
        entry = ttk.Entry(masterFrame, textvariable=variable, width=10)
    entry.grid(row=row, column=column, sticky="W")
    variable.set(default)

    labelUnit = None
    if (unit != ""):
        labelUnit = ttk.Label(masterFrame, text="["+unit+"]")
        labelUnit.grid(row=row, column=0, sticky="E")
    return(label, entry, labelUnit)


#------------------------------------------------------------------------------*
# Widget Separator

def seperator(frame, row: int, column: int=0, orientation: str="h") -> int:
    """
    Widget seperator line. Can be horizontal using 'h' or vertical using 'v' in
    the parameter 'orientation'.
    """
    if (orientation == "h"):
        ttk.Separator(frame, orient=tk.HORIZONTAL).grid(column=column, columnspan=10, row=row, sticky='we')
    elif (orientation == "v"):
        ttk.Separator(frame, orient=tk.VERTICAL).grid(column=column, row=row, rowspan=100, sticky='ns')
    return(row+1)


#------------------------------------------------------------------------------*
# Dropdown selector from a list

def dropdownSelection(
    frame,
    selectionVariable,
    textLabel: str,
    variableList: dict,
    row: int,
    column: int=2
) -> int:
    ttk.Label(frame, text=textLabel).grid(row=row, column=0, sticky="W")

    # Set the default option
    selectionVariable.set(variableList[0])
    variableList.insert(0, "")
    dropdown = ttk.OptionMenu(frame, selectionVariable, *variableList)
    dropdown.grid(row=row, column=column, sticky="W")

    return(row+1)

#------------------------------------------------------------------------------*
# Image

def logo(frame):
    byte_data = base64.b64decode(image.imageByte)
    image_data = io.BytesIO(byte_data)
    img = Image.open(image_data)

    # Label
    img = ImageTk.PhotoImage(Image.open(image_data).resize((80, 63), Image.ANTIALIAS))
    labelLogo = ttk.Label(frame, image=img)
    labelLogo.image = img

    return(labelLogo, img)
