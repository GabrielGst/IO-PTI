"""
This module provides functions for creating and customizing plots using matplotlib. Default export PATH is set to current directory.

Functions:
    saveImsubplots(data, legends, filename="output", columnNumber=3, size=(12, 4), save=False, plot=True):
        Save and plot subplots of images or arrays with plt.imshow()
"""

from math import ceil
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

SAVE_PATH = "./"

def savePlot(x, y, title, filename="output", save=False, plot=True, logX=False, logY=False):
    """Save and plot a single serie with plt.plot()

    Args:
        x (array): x-axis values
        y (array): y-axis values
        title (string): Title of the plot
        filename (string, optional): Filename to save the figure. Defaults to "output".
        save (bool, optional): Defaults to False.
        plot (bool, optional): Defaults to True.
    """
    
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_title(title)
    ax.grid(True, linestyle='-.')

    if logX:
        ax.set_xscale('log')
    if logY:
        ax.set_yscale('log')
    
    if save:
        fig.savefig(f'{SAVE_PATH}' + f'{filename}.png') # Save the figure for remote shell
    if plot:
        plt.show()
    else:
        print("Please select save or plot option")


def saveImplot(data, title, filename="output", save=False, plot=True, logScale=False):
    """Save and plot a single image or array with plt.imshow()

    Args:
        data (array): Array to plot
        title (string): Title of the plot
        filename (string, optional): Filename to save the figure. Defaults to "output".
        save (bool, optional): Defaults to False.
        plot (bool, optional): Defaults to True.
    """
    
    plt.figure()
    plt.imshow(data, cmap='gray', norm=LogNorm() if logScale else None)
    plt.title(title)
    plt.colorbar()
    if save:
        plt.savefig(f'{SAVE_PATH}' + f'{filename}.png') # Save the figure for remote shell
    if plot:
        plt.show()
    else:
        print("Please select save or plot option")

def saveImsubplots(data, legends, filename="output", columnNumber=3, size=(12, 4), save=False, plot=True, logScale=False):
    """Save and plot subplots of images or arrays with plt.imshow()

    Args:
        data (list): List of arrays to plot
        legends (list): List of legends for each subplot
        filename (string, optional): Filename to save the figure. Defaults to "output".
        columnNumber (int, optional): Number of co:umns. Defaults to 3.
        figsize (tuple, optional): Tuple of 2 ints, (width, height). Defaults to (12, 4).
        save (bool, optional): Defaults to False.
        plot (bool, optional): Defaults to True.
    """
    
    n = len(legends)
    rowNumber = ceil(n / columnNumber)
    plt.figure(figsize=size) # Create a 1x3 grid for the subplots

    for i in range(n):
        # n-th subplot
        plt.subplot(rowNumber, columnNumber, i+1)  # (rows, columns, index)
        plt.imshow(data[i], cmap='gray', norm=LogNorm() if logScale else None)
        plt.colorbar()
        plt.title(f'{legends[i]}')

    # Adjust layout for better spacing
    plt.tight_layout()

    # Show the plots
    if save:
        plt.savefig(f'{SAVE_PATH}' + f'{filename}.png') # Save the figure for remote shell
    if plot:
        plt.show()
    else:
        print("Please select save or plot option")
