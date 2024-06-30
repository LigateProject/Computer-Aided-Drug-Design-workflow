#!/usr/bin/python3

"""
This is a primitive replacement for uncertainties.ufloat. It covers only
the most basic functionality.

Author: Magnus Lundborg
License: Apache License, version 2.0
"""

import math

class internal_ufloat:
    """
    The class is used to propagate uncertainties of float values
    Currently only used for simple addition and subtraction.
    No error correlation is taken into account.

    The uncertainty is called std_dev to be compliant with
    uncertainties.ufloat.

    For full functionality use the uncertainties.ufloat class instead.
    """

    def __init__(self, value, std_dev = None):
        if isinstance(value, internal_ufloat):
            self.nominal_value = value.nominal_value
            if std_dev is not None:
                self.std_dev = std_dev
            else:
                self.std_dev = value.std_dev
        elif isinstance(value, str):
            parts = value.split('+/-')
            self.nominal_value = float(parts[0])
            self.std_dev = float(parts[1])
        else:
            self.nominal_value = float(value)
            self.std_dev = std_dev

    def __repr__(self):
        if self.std_dev is not None:
            return f'{self.nominal_value:.2f}+/-{self.std_dev:.2f}'
        return f'{self.nominal_value:.2f}'

    def __str__(self):
        if self.std_dev is not None:
            return f'{self.nominal_value:.2f}+/-{self.std_dev:.2f}'
        return f'{self.nominal_value:.2f}'

    def __add__(self, other):
        if isinstance(other, internal_ufloat):
            std_dev = math.sqrt(self.std_dev ** 2 + other.std_dev ** 2)
            return internal_ufloat(self.nominal_value + other.nominal_value, std_dev)

        std_dev = self.std_dev
        return internal_ufloat(self.nominal_value + float(other), std_dev)

    def __radd__(self, other):
        if isinstance(other, internal_ufloat):
            std_dev = math.sqrt(self.std_dev ** 2 + other.std_dev ** 2)
            return internal_ufloat(other.nominal_value + self.nominal_value, std_dev)

        return internal_ufloat(float(other) + self.nominal_value, self.std_dev)

    def __sub__(self, other):
        if isinstance(other, internal_ufloat):
            std_dev = math.sqrt(self.std_dev ** 2 + other.std_dev ** 2)
            return internal_ufloat(self.nominal_value - other.nominal_value, std_dev)

        return internal_ufloat(self.nominal_value - float(other), self.std_dev)

    def __rsub__(self, other):
        if isinstance(other, internal_ufloat):
            std_dev = math.sqrt(self.std_dev ** 2 + other.std_dev ** 2)
            return internal_ufloat(other.nominal_value - self.nominal_value, std_dev)

        return internal_ufloat(float(other) - self.nominal_value, self.std_dev)

    def __lt__(self, other):
        if isinstance(other, internal_ufloat):
            return self.nominal_value < other.nominal_value

        return self.nominal_value < float(other)

    def __gt__(self, other):
        if isinstance(other, internal_ufloat):
            return self.nominal_value > other.nominal_value

        return self.nominal_value > float(other)

    def __iadd__(self, other):
        if isinstance(other, internal_ufloat):
            self.nominal_value = self.nominal_value + other.nominal_value
            self.std_dev = math.sqrt(self.std_dev ** 2 + other.std_dev ** 2)
        else:
            self.nominal_value = self.nominal_value + float(other)
        return internal_ufloat(self.nominal_value, self.std_dev)

    def __isub__(self, other):
        if isinstance(other, internal_ufloat):
            self.nominal_value = self.nominal_value - other.nominal_value
            self.std_dev = math.sqrt(self.std_dev ** 2 + other.std_dev ** 2)
        else:
            self.nominal_value = self.nominal_value - float(other)
        return internal_ufloat(self.nominal_value, self.std_dev)
