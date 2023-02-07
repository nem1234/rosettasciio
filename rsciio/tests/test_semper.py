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

import os.path
import tempfile
from os import remove

import numpy as np
import pytest

hs = pytest.importorskip("hyperspy.api", reason="hyperspy not installed")

my_path = os.path.dirname(__file__)

# Reference data:
data_signal = np.arange(27, dtype=np.float32).reshape((3, 3, 3)) / 2.0
data_image = np.arange(16, dtype=np.float32).reshape((4, 4)) / 2.0
data_spectrum = np.arange(10, dtype=np.float32) / 2.0
data_image_byte = np.arange(25, dtype=np.byte).reshape(
    (5, 5)
)  # Odd dim. tests strange read/write
data_image_int16 = np.arange(16, dtype=np.int16).reshape((4, 4))
data_image_int32 = np.arange(16, dtype=np.int32).reshape((4, 4))
data_image_complex = (data_image_int32 + 1j * data_image).astype(np.complex64)
test_title = "This is a test!"


def test_writing_unsupported_data_type():
    data = np.arange(5 * 10).reshape((5, 10))
    s = hs.signals.BaseSignal(data.astype("int64"))
    with tempfile.TemporaryDirectory() as tmpdir:
        with pytest.raises(IOError) as cm:
            fname = os.path.join(tmpdir, "test_writing_unsupported_data_type.unf")
            s.save(fname)
            cm.match(
                "The SEMPER file format does not support int64 data type",
            )


def test_writing_loading_metadata():
    data = np.arange(5 * 10).reshape((5, 10)).astype(np.int8)
    s = hs.signals.BaseSignal(data)
    s.metadata.set_item("General.date", "2016-08-06")
    s.metadata.set_item("General.time", "11:55:00")
    with tempfile.TemporaryDirectory() as tmpdir:
        fname = os.path.join(tmpdir, "test_write_with_metadata.unf")
        s.save(fname)
        s2 = hs.load(fname)
        np.testing.assert_allclose(s.data, s2.data)
        assert s.metadata.General.date == s2.metadata.General.date
        assert s.metadata.General.time == s2.metadata.General.time


def test_signal_3d_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_signal_3d.unf"))
    np.testing.assert_equal(signal.data, data_signal)
    np.testing.assert_equal(signal.original_metadata.IFORM, 2)  # float
    assert isinstance(signal, hs.signals.BaseSignal)


def test_image_2d_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_image_2d.unf"))
    np.testing.assert_equal(signal.data, data_image)
    np.testing.assert_equal(signal.original_metadata.IFORM, 2)  # float
    assert isinstance(signal, hs.signals.Signal2D)


def test_spectrum_1d_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_spectrum_1d.unf"))
    np.testing.assert_equal(signal.data, data_spectrum)
    np.testing.assert_equal(signal.original_metadata.IFORM, 2)  # float
    assert isinstance(signal, hs.signals.Signal1D)


def test_image_byte_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_image_byte.unf"))
    np.testing.assert_equal(signal.data, data_image_byte)
    np.testing.assert_equal(signal.original_metadata.IFORM, 0)  # byte
    assert isinstance(signal, hs.signals.Signal2D)


def test_image_int16_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_image_int16.unf"))
    np.testing.assert_equal(signal.data, data_image_int16)
    np.testing.assert_equal(signal.original_metadata.IFORM, 1)  # int16
    assert isinstance(signal, hs.signals.Signal2D)


def test_image_int32_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_image_int32.unf"))
    np.testing.assert_equal(signal.data, data_image_int32)
    np.testing.assert_equal(signal.original_metadata.IFORM, 4)  # int32
    assert isinstance(signal, hs.signals.Signal2D)


def test_image_complex_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_image_complex.unf"))
    np.testing.assert_equal(signal.data, data_image_complex)
    np.testing.assert_equal(signal.original_metadata.IFORM, 3)  # complex
    assert isinstance(signal, hs.signals.ComplexSignal)


def test_with_title_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_with_title.unf"))
    np.testing.assert_equal(signal.data, data_image)
    np.testing.assert_equal(signal.original_metadata.IFORM, 2)  # float
    np.testing.assert_equal(signal.metadata.General.title, test_title)
    assert isinstance(signal, hs.signals.Signal2D)


def test_no_label_loading():
    signal = hs.load(os.path.join(my_path, "semper_data", "example_no_label.unf"))
    np.testing.assert_equal(signal.data, data_image)
    np.testing.assert_equal(signal.original_metadata.ILABEL, 0)
    assert isinstance(signal, hs.signals.Signal2D)


class TestCaseSaveAndReadImage:
    def test_save_and_read(self):
        signal_ref = hs.signals.Signal2D(data_image)
        signal_ref.metadata.General.title = test_title
        signal_ref.save(
            os.path.join(my_path, "semper_data", "example_temp.unf"), overwrite=True
        )
        signal = hs.load(os.path.join(my_path, "semper_data", "example_temp.unf"))
        np.testing.assert_equal(signal.data, signal_ref.data)
        np.testing.assert_equal(signal.metadata.General.title, test_title)
        assert isinstance(signal, hs.signals.Signal2D)

    def teardown_method(self, method):
        remove(os.path.join(my_path, "semper_data", "example_temp.unf"))


class TestCaseSaveAndReadByte:
    def test_save_and_read(self):
        signal_ref = hs.signals.Signal2D(data_image_byte)
        signal_ref.metadata.General.title = test_title
        signal_ref.save(
            os.path.join(my_path, "semper_data", "example_temp.unf"), overwrite=True
        )
        signal = hs.load(os.path.join(my_path, "semper_data", "example_temp.unf"))
        np.testing.assert_equal(signal.data, signal_ref.data)
        np.testing.assert_equal(signal.metadata.General.title, test_title)
        assert isinstance(signal, hs.signals.Signal2D)

    def teardown_method(self, method):
        remove(os.path.join(my_path, "semper_data", "example_temp.unf"))
