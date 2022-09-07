# %%
# -*- coding: utf-8 -*-
#
# helpers.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

r"""Helper functions for the Sudoku solver
----------------------------------------------------------------

:Authors: J Gille, S Furber, A Rowley
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch


# RGB values for the colors used in the generated images
green = (0, 115, 0)
red = (180, 0, 0)
black = (0, 0, 0)
dark_grey = (50, 50, 50)
light_grey = (160, 160, 160)
white = (255, 255, 255)


def get_puzzle(puzzle_index):
    """returns one of 8 Sudoku configuration to be solved.

    Args:
        puzzle_index (int): index between 0 and 7 indicating the puzzle number

    Returns:
        np.array: array of shape (9,9) representing the puzzle configuration.
        Array is zero wherever no input is given, and contains the corresponding
        digit otherwise.
    """
    init_config = None

    if not 0 <= puzzle_index < 8:
        raise ValueError(
            "Cannot return puzzle - index must be between 0 and 7!")

    if puzzle_index == 0:
        # Dream problem: make the network come up with a valid sudoku without
        # any restrictions
        init_config = np.zeros((9, 9), dtype=np.uint8)
    elif puzzle_index == 1:
        # Diabolical problem:
        init_config = [[0, 0, 1,  0, 0, 8,  0, 7, 3],
                       [0, 0, 5,  6, 0, 0,  0, 0, 1],
                       [7, 0, 0,  0, 0, 1,  0, 0, 0],

                       [0, 9, 0,  8, 1, 0,  0, 0, 0],
                       [5, 3, 0,  0, 0, 0,  0, 4, 6],
                       [0, 0, 0,  0, 6, 5,  0, 3, 0],

                       [0, 0, 0,  1, 0, 0,  0, 0, 4],
                       [8, 0, 0,  0, 0, 9,  3, 0, 0],
                       [9, 4, 0,  5, 0, 0,  7, 0, 0]]

    elif puzzle_index == 2:
        init_config = [[2, 0, 0,  0, 0, 6,  0, 3, 0],
                       [4, 8, 0,  0, 1, 9,  0, 0, 0],
                       [0, 0, 7,  0, 2, 0,  9, 0, 0],

                       [0, 0, 0,  3, 0, 0,  0, 9, 0],
                       [7, 0, 8,  0, 0, 0,  1, 0, 5],
                       [0, 4, 0,  0, 0, 7,  0, 0, 0],

                       [0, 0, 4,  0, 9, 0,  6, 0, 0],
                       [0, 0, 0,  6, 4, 0,  0, 1, 9],
                       [0, 5, 0,  1, 0, 0,  0, 0, 8]]

    elif puzzle_index == 3:
        init_config = [[0, 0, 3,  2, 0, 0,  0, 7, 0],
                       [0, 0, 5,  0, 0, 0,  3, 0, 0],
                       [0, 0, 8,  9, 7, 0,  0, 5, 0],

                       [0, 0, 0,  8, 9, 0,  0, 0, 0],
                       [0, 5, 0,  0, 0, 0,  0, 2, 0],
                       [0, 0, 0,  0, 6, 1,  0, 0, 0],

                       [0, 1, 0,  0, 2, 5,  6, 0, 0],
                       [0, 0, 4,  0, 0, 0,  8, 0, 0],
                       [0, 9, 0,  0, 0, 7,  5, 0, 0]]

    elif puzzle_index == 4:
        init_config = [[0, 1, 0,  0, 0, 0,  0, 0, 2],
                       [8, 7, 0,  0, 0, 0,  5, 0, 4],
                       [5, 0, 2,  0, 0, 0,  0, 9, 0],

                       [0, 5, 0,  4, 0, 9,  0, 0, 1],
                       [0, 0, 0,  7, 3, 2,  0, 0, 0],
                       [9, 0, 0,  5, 0, 1,  0, 4, 0],

                       [0, 2, 0,  0, 0, 0,  4, 0, 8],
                       [4, 0, 6,  0, 0, 0,  0, 1, 3],
                       [1, 0, 0,  0, 0, 0,  0, 2, 0]]

    elif puzzle_index == 5:
        init_config = [[8, 9, 0,  2, 0, 0,  0, 7, 0],
                       [0, 0, 0,  0, 8, 0,  0, 0, 0],
                       [0, 4, 1,  0, 3, 0,  5, 0, 0],

                       [2, 5, 8,  0, 0, 0,  0, 0, 6],
                       [0, 0, 0,  0, 0, 0,  0, 0, 0],
                       [6, 0, 0,  0, 0, 0,  1, 4, 7],

                       [0, 0, 7,  0, 1, 0,  4, 3, 0],
                       [0, 0, 0,  0, 2, 0,  0, 0, 0],
                       [0, 2, 0,  0, 0, 7,  0, 5, 1]]

    elif puzzle_index == 6:
        # "World's hardest sudoku":
        # http://www.telegraph.co.uk/news/science/science-news/9359579/Worlds-hardest-sudoku-can-you-crack-it.html
        init_config = [[8, 0, 0,  0, 0, 0,  0, 0, 0],
                       [0, 0, 3,  6, 0, 0,  0, 0, 0],
                       [0, 7, 0,  0, 9, 0,  2, 0, 0],

                       [0, 5, 0,  0, 0, 7,  0, 0, 0],
                       [0, 0, 0,  0, 4, 5,  7, 0, 0],
                       [0, 0, 0,  1, 0, 0,  0, 3, 0],

                       [0, 0, 1,  0, 0, 0,  0, 6, 8],
                       [0, 0, 8,  5, 0, 0,  0, 1, 0],
                       [0, 9, 0,  0, 0, 0,  4, 0, 0]]

    elif puzzle_index == 7:
        init_config = [[1, 0, 0,  4, 0, 0,  0, 0, 0],
                       [7, 0, 0,  5, 0, 0,  6, 0, 3],
                       [0, 0, 0,  0, 3, 0,  4, 2, 0],

                       [0, 0, 9,  0, 0, 0,  0, 3, 5],
                       [0, 0, 0,  3, 0, 5,  0, 0, 0],
                       [6, 3, 0,  0, 0, 0,  1, 0, 0],

                       [0, 2, 6,  0, 5, 0,  0, 0, 0],
                       [9, 0, 4,  0, 0, 6,  0, 0, 7],
                       [0, 0, 0,  0, 0, 8,  0, 0, 2]]

    return np.array(init_config)


def validate_solution(puzzle, solution):
    """validate a proposed solution for a sudoku puzzle

    Args:
        puzzle (np.array): array of shape (9,9) encoding the puzzle.
        see get_puzzle().
        solution (np.array): array of shape (9,9) encoding the proposed
        solution.

    Returns:
        (bool, np.array, np.array, np.array): tuple of values that indicate
        the validity of the solution:
        1. True if the overall solution is valid, False otherwise.
        2. boolean array of shape (3,3) that is True wherever a 3x3 box is valid
        3. boolean array of shape (9,) encoding the validity of all rows
        4. boolean array of shape (9,) encoding the validity of all columns
    """

    boxes = np.ones((3, 3), dtype=bool)
    rows = np.ones(9, dtype=bool)
    cols = np.ones(9, dtype=bool)

    expected_numbers = set(range(1, 10))

    # validate boxes
    for i in range(3):
        for j in range(3):
            box = solution[3*i:3*i+3, 3*j:3*j+3]
            if expected_numbers != set(box.flatten()):
                boxes[i, j] = False

    # validate rows and columns
    for i in range(9):
        if expected_numbers != set(solution[i, :]):
            rows[i] = False
        if expected_numbers != set(solution[:, i]):
            cols[i] = False

    # It is possible (in rare cases) that the network finds a valid
    # solution that does not conform to the initial puzzle configuration, that is,
    # one of the cells where input is applied is overridden by the rest of the
    # network. This is taken care of here.
    input_cells = np.where(puzzle != 0)
    puzzle_matched = puzzle[input_cells] == solution[input_cells]

    # overall solution is valid iff all of the components are.
    valid = boxes.all() and rows.all() and cols.all() and puzzle_matched.all()

    return valid, boxes, rows, cols


cell_size = 18  # inner size of the cell
grid_width = 2  # width of the grid separating cells
# step size from a position in a cell to the same position in the next
cell_step = cell_size + grid_width

field_size = 9 * cell_size + 10 * grid_width
frame_width = 10
image_size = field_size + 2*frame_width


def plot_field(puzzle, solution, with_color=False):
    """generates a graphical representation of a Sudoku field. Digits that are
    given by the puzzle are represented as bold and black, while calculated
    digits are represented in grey and italic.

    Args:
        puzzle (np.array): array of shape (9,9) that represents the puzzle
        that is being solved. See get_puzzle()
        solution (np.array): array of shape (9,9) representing the solution.
        with_color (bool, optional): if True, green and red are used to
        indicate which parts of the solution are valid and which are not.
        Otherwise, only black and white are used. Defaults to False.

    Returns:
        PIL.Image: A visual representation of the Sudoku solution.
    """
    ax2 = plt.subplot2grid((20, 3), (10, 2), rowspan=10, colspan=1)
    decorate_sudoku_box(ax2)
    fill_numbers(ax2, puzzle, solution)

    # color every other box in light grey
    for i in range(0, 9, 3):
        for j in range(0, 9, 3):
            if (i+j) % 2 == 0:
                ax2.add_patch(patch.Rectangle((j, i), 3, 3, facecolor="grey",
                                              alpha=0.5))
    if with_color:
        valid, boxes, rows, cols = validate_solution(puzzle, solution)
        # draw a red or green line in the background to indicate validity of
        # all rows and columns

        for i in range(9):
            ax2.vlines(np.where(rows)[0] +0.5, -0.5, 9.5, color="red")
            ax2.hlines(np.where(cols), 0, 9, color="red")

        # draw a red frame around boxes that are not valid
        for i in range(3):
            for j in range(3):
                if not boxes[i, j]:
                    print(f"patching: {i,j}")
                    ax2.add_patch(patch.Rectangle((3*j, 3*j + 0.1), 2.9, 2.9,
                                                  color="red", fill=False, linewidth=3))
    return ax2


def fill_numbers(axis, puzzle, solution):
    # fill the numbers
    for i in range(9):
        for j in range(9):
            # quirk: when plotting matrix, transpose it; i is in y-axis and j is in x-axis
            if puzzle[i][j] == 0:
                text_style = 'normal'
                text_col = 'black'
            else:
                text_style = 'italic'
                if puzzle[i][j] == solution[i][j]:
                    text_col = 'gray'
                else:
                    # If the network proposes a solution where a digit
                    # from the input configuration is altered, that
                    # digit is colored in red.
                    text_col = 'red'

            axis.text(j + 0.5, i + 0.5, solution[i][j], horizontalalignment='center',
                      verticalalignment='center', color=text_col, fontsize=9, style=text_style)


def decorate_sudoku_box(axis):
    """
    Modify line and color properties of sudoku grid to make it look good
    """
    [x, y, xr, yr] = generate_sudoku_box_lines()
    axis.plot(x, y, color='gray')
    axis.plot(xr, yr, color='gray', linewidth=2)
    # turn off tick markers
    axis.set_xticks([])
    axis.set_yticks([])



def generate_sudoku_box_lines():
    # TODO: unfuck this shite!
    # lines for cells
    x = []
    y = []
    # lines for regions
    xr = []
    yr = []

    # data for vertical lines
    for i in range(10):
        x.append(i)
        x.append(i)
        y.append(0)
        y.append(9)
        x.append(None)
        y.append(None)
        if divmod(i, 3)[1] == 0:
            xr.append(i)
            xr.append(i)
            yr.append(0)
            yr.append(9)
            xr.append(None)
            yr.append(None)

    # data for horizontal lines
    for j in range(10):
        x.append(0)
        x.append(9)
        y.append(j)
        y.append(j)
        x.append(None)
        y.append(None)
        if divmod(j, 3)[1] == 0:
            xr.append(0)
            xr.append(9)
            yr.append(j)
            yr.append(j)
            xr.append(None)
            yr.append(None)


    return [x, y, xr, yr]


# %%
