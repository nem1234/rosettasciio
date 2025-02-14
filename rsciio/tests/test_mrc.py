# -*- coding: utf-8 -*-
# Copyright 2007-2023 The HyperSpy developers
#
# This file is part of RosettaSciIO.
#
# RosettaSciIO is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RosettaSciIO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with RosettaSciIO. If not, see <https://www.gnu.org/licenses/#GPL>.

from pathlib import Path
import pytest

hs = pytest.importorskip("hyperspy.api", reason="hyperspy not installed")


TEST_DATA_DIR = Path(__file__).parent / "data" / "mrc"


def test_single_image():
    # Acquired from Velox
    s = hs.load(TEST_DATA_DIR / "HAADFscan.mrc")
    assert s.data.shape == (16, 16)
    assert s.axes_manager.signal_shape == (16, 16)
    assert s.axes_manager.navigation_shape == ()

    for axis in s.axes_manager.signal_axes:
        assert axis.scale == 5.679131317138672
        assert axis.offset == 0
        assert axis.units == "nm"


def test_4DSTEM_image():
    # Acquired from Velox
    s = hs.load(TEST_DATA_DIR / "4DSTEMscan.mrc")
    assert s.data.shape == (256, 256, 256)
    assert s.axes_manager.signal_shape == (256, 256)
    assert s.axes_manager.navigation_shape == (256,)


def test_4DSTEM_image_navigation_shape_16_16():
    # Acquired from Velox
    s = hs.load(
        TEST_DATA_DIR / "4DSTEMscan.mrc",
        navigation_shape=(16, 16),
    )
    assert s.data.shape == (16, 16, 256, 256)
    assert s.axes_manager.signal_shape == (256, 256)
    assert s.axes_manager.navigation_shape == (16, 16)


def test_4DSTEM_image_navigation_shape_8_32():
    s = hs.load(
        TEST_DATA_DIR / "4DSTEMscan.mrc",
        navigation_shape=(8, 32),
    )
    assert s.data.shape == (32, 8, 256, 256)
    assert s.axes_manager.signal_shape == (256, 256)
    assert s.axes_manager.navigation_shape == (8, 32)
