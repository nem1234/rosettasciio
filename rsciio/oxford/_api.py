# -*- coding: utf-8 -*-
#
#  ================== DRAFT ==============
#  Not ready to merge into main branch
#  File reader for Oxford EDS/EBSD system
#
# Supported:
#   read SEM Image with metadata
#   read EDS Map with metadata
# ToDo:
#   read EDS Spectrum
#   read EBSD Image
#   reduce mem size
#   lazy loading
#
#
# Copyright 2023- The HyperSpy developers
#
# This file will be part of RosettaSciIO.
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

import os
from datetime import datetime, timedelta, timezone
import logging

import copy
import struct
import sys

import sqlite3
from pathlib import Path
# from box import Box, BoxList

import numpy as np
from rsciio.utils import tools
from rsciio._docstrings import FILENAME_DOC, LAZY_DOC, RETURNS_DOC



_logger = logging.getLogger(__name__)

class Utils:
    def _get_simple_binary_block(self, length):
        _data = self.data[0:length]
        self._rseek(length)
        return _data

    def _get_string(self, data=None):
        if data is None:
            byte_length = self.data[0]
            _data = self.data[1:1 + byte_length]
            self._rseek(1 + byte_length)
            return _data.decode()
        else:
            byte_length = data[0]
            return data[1:1 + byte_length], data[1 + byte_length:]



    def _get_num(self, type_str, data):
        _len = struct.calcsize(type_str)
        if data is None:
            _data = struct.unpack(type_str, self.data[0:_len])
            self._rseek(_len)
        else:
            _data = struct.unpack(type_str, data[0:_len])
            data = data[_len:]
        if data is None:
            if len(_data) == 1:
                return _data[0]
            return _data
        else:
            if len(_data) == 1:
                return _data[0], data
            return _data, data

    def _get_type(self, data=None):
        return self._get_num('<B', data)

    def _get_byte(self, data=None):
        return self._get_num('<B', data)

    def _get_ulong(self, data=None):
        return self._get_num('<L', data)

    def _get_long(self, data=None):
        return self._get_num('<l', data)

    def _get_ushort(self, data=None):
        return self._get_num('<H', data)

    def _get_short(self, data=None):
        return self._get_num('<h', data)

    def _get_float(self, data=None):
        return self._get_num('<f', data)

    def _get_double(self, data=None):
        return self._get_num('<d', data)

    def __init__(self):
        self.data = b''

    def _dump_bytes_debug(self, data=None):
        if data is None:
            _len = len(self.data)
            if _len > 16:
                _len = 16
            data = self.data[0:16]
        string = ""
        for ch in data:
            string += f'{np.uint8(ch):02x} '
        return string


class _BinaryMetadata(Utils):
    def __init__(self):
        """Parameters:
        ----------------
        data : bytes
            internal bynary data
        """
        self.original_data = None
        self.data = None
        self.datanodes = None
        self.debug = 2

    def _tell(self):
        return len(self.original_data)-len(self.data)

    def _seek(self, offset):
        if offset < 0 or offset > len(self.original_data):
            print("Seek offset = ", offset, '/', len(self.original_data))
            _logger.warning('Data length is shorter than offset')
        self.data = self.original_data[offset:]
    def _rseek(self, offset):
#        print("rseek offset = ", offset)
        self._seek(self._tell() + offset)

    def _decode_metadata_internal(self):
        """
        Read binary data (serialized XML?)

        Parameters
        ----------------------------
        Returns:
                metadata : dict
        """
        block = { "name": "", "children": [], "values": []}
        namespace = ""
        scheme = ""
        block_type = None
        subblock_name = None
        # for debug
        with open("test.bin", "wb") as f:
            f.write(self.data)

        block_type = self._get_type()
#        if block_type != 0x40:
        if block_type not in (0x40, 0x5e, 0x5f, 0x61):
            print(f"Unknown type: {block_type:02x}")
            return {}
        block["name"] = self._get_string()
        print(self.data)
        while len(self.data) > 0:
            ptr = self._tell() # keep original pointer
            block_type = self._get_type()
            if block_type == 0x01:
                if subblock_name is not None:
                    subblock_name = None
                else:
                    break
            elif block_type in (0x40, 0x5e, 0x5f, 0x61): # Block
                block_name = self._get_string()
                if block_name == "Blob":
                    block["values"].append(self._get_binary_block())
                elif self.data[0] in (0x8b, 0x97):
                    print(block_name)
                    _t = self._get_binary_block()
                    block[block_name]=_t
                else:
                    self._seek(ptr)
                    print("child=", self._dump_bytes_debug())
                    _child = self._decode_metadata_internal()
                    block["children"].append(_child)
            elif block_type in (0x5e, 0x5f, 0x61):
                subblock_name = self._get_string()
                ptr = self._tell()
                block_type = self._get_type()
                #if block_type == 0x2e:
                #    # skip data type
                #    _attr_type = self._get_string()
                #    _attr_data = self._get_binary_block()
                _value = self._get_binary_block()
                if subblock_name not in block:
                    block[subblock_name] = []
                if not isinstance(block[subblock_name], list):
                    block[subblock_name] = [block[subblock_name]]
            elif block_type == 0x60:
                block_name = "struct"
                __data_type = self._get_string()  # "anyData"
            elif block_type == 0x2e:  # Attribute
                _typ = self._get_string()
                if _typ == "type":
                    self._rseek(1)
                    _type = self._get_string()
                    block[_typ] = _type
            elif block_type == 0x08:   # NameSpace ?
                namespace = self._get_string()
                block["namespace"] = namespace
            elif block_type == 0x09: # Scheme
                _typ_unknown = self._get_type()
                _typ2 =  self._get_type()
                scheme = self._get_string()
                block["scheme"]  = scheme
            elif block_type == 0x86: # unknown, skip 1 byte
                self._get_byte()
            else:
                print(f"Unknown {block_type:02x} ", self._dump_bytes_debug())
        return block

    def _get_binary_block(self, block_type=None):
        if block_type is None:
            block_type = self._get_type()
        _blk = b''
        while True:
            if block_type in (0x9e, 0x9f,0x98):
                _blen = self._get_type()
                _blk += self._get_simple_binary_block(_blen)
            elif block_type == 0xa0:
                _blen = self._get_ushort()
                _blk += self._get_simple_binary_block(_blen)
            elif block_type == 0x97: # Date Time 8 bytes ??
                _blk = self._get_simple_binary_block(8)
            elif block_type == 0x8b: # OffsetMinutes (TZ), 3 bytes
                _blk = self._get_short()
            elif block_type == 0x86: # unknown, skip 1 bytes
                _ = self.get_byte()
            if block_type in (0x86, 0x8b, 0x97, 0x98, 0x9f):  # Single block data or End of data
                break
            if block_type == 1:  # End of data block?
                break
            block_type = self._get_type()
        if block_type == 0x98: # convert to text
            _blk = _blk.decode()
        return _blk



    def decode_metadata(self, row, columns):
        """ convert from binary data in db recored to text metadata
        """
        print(columns)
        self.original_data = row[columns["Data"]]
        self.data = self.original_data
        _dict = self._decode_metadata_internal()
        return _dict




class MetaData(Utils):
    def _convert_metadata_to_dict_internal(self, dic):
        for _key, _item in list(dic.items()):
            if _key == "Blob":
                dic["values"].append(_item)
            elif _key == "children":
                for _child in _item:
                    _childdict = self._convert_metadata_to_dict_internal(_child)
                    dic[_child["name"]] = _childdict
        if "Blob" in dic:
            del dic["Blob"]
        del dic["children"]
        return dic

    def _convert_metadata_to_dict(self, raw_dict):
        _dict = self._convert_metadata_to_dict_internal(raw_dict)
        return _dict


    def _decode_uuid(self, bin_uuid):
        _b = bin_uuid
        _uuid = f'{_b[3]:02x}{_b[2]:02x}{_b[1]:02x}{_b[0]:02x}-' \
            f'{_b[5]:02x}{_b[4]:02x}-{_b[7]:02x}{_b[6]:02x}-' \
            f'{_b[8]:02x}{_b[9]:02x}-{_b[10]:02x}{_b[11]:02x}' \
            f'{_b[12]:02x}{_b[13]:02x}{_b[14]:02x}{_b[15]:02x}'
        return _uuid

    def _get_dict(self, data, dtype='<d', dlen=8):
        # decode parameter list paired with its name
        _n_data, data = self._get_long(data)
        _keys = []
        _vals = []
        for i in range(_n_data):
            _len, data = self._get_byte(data[1:])
            _keys.append(data[0:_len].decode())
            data = data[_len:]
        _, data = self._get_long(data)  # should be same as _n_data
        for i in range(_n_data):
            _val = struct.unpack(dtype, data[0:dlen])[0]
            data = data[dlen:]
            _vals.append(_val)
        return dict(zip(_keys, _vals)), data

    #
    # record type specific decoding functions
    #
    def _post_proc_Project(self, dic, **kwargs):
        # Currently only TimeZone data is used.
        new_dic = {}
        new_dic['UUID'] = self._decode_uuid(dic["values"][0][4:20])
        new_dic['DataCounters'], _ = self._get_dict(dic['DataCounters']['values'][0][4:], '<l', 4)
        new_dic['OffsetMinutes'] = dic['CreatedOn']['OffsetMinutes']
        self.time_zone = dic['CreatedOn']['OffsetMinutes']
        dic = self._post_proc_cleanup(new_dic)
        return dic

    def _post_proc_DataNode(self, dic):
        _dat = dic["values"][0]
        new_dic = tools.DTBox(box_dots=True, default_box=True)
        uuid = self._decode_uuid(_dat[4:20])
        new_dic['UUID'] = uuid
        if _dat[24] == 0:
            new_dic['Text'], _dat = self._get_string(_dat[25:])
        else:
            _dat = _dat[25:]
        new_dic['UUID_parent'] = self._decode_uuid(_dat[0:16])
        new_dic['UUID_data'] = self._decode_uuid(_dat[16:32])
        new_dic['Type'] = struct.unpack('<l', _dat[32:36])[0]
        new_dic['Unknown1'] = _dat[36]
        new_dic['SortKey'] = struct.unpack('<l', _dat[37:41])[0]
        _have_text2 = struct.unpack('<l', _dat[41:45])[0]
        if _have_text2:
            new_dic["Text2"], _dat = self._get_string(_dat[46:])
        else:
            _dat = _dat[45:]
        if len(_dat) > 0:
            new_dic['Unknown3'] = self._dump_bytes_debug(_dat)
        return {uuid: new_dic}


    def _post_proc_cleanup(self, dic, dellist=[]):
        _dellist = []
        if isinstance(dic, dict):
            for _key, _val in dic.items():
                if _key in dellist:
                    _dellist.append(_key)
                if isinstance(_val, list) or isinstance(_val, dict):
                    self._post_proc_cleanup(_val, dellist)
                if _key == 'values' and isinstance(_val, list):
                    _dellist.append(_key)
                    # dic[_key] = [self._dump_bytes_debug(_v[0:16]) for _v in _val]
            _dellist.extend(dellist)
            print(_dellist)
            for _itm in _dellist:
                print("Deleting ", _itm)
                if _itm in dic:
                    del dic[_itm]
        elif isinstance(dic, list):
            for _val in dic:
                if isinstance(_val, list) or isinstance(_val, dict):
                    self._post_proc_cleanup(_val, dellist)
        return dic


    def _post_proc_AreaOfInterest(self, dic, **kwargs):
        new_dic = tools.DTBox(box_dots=True, default_box=True)
        _uuid = self._decode_uuid(dic["values"][0][4:20])
        new_dic["UUID"] = _uuid
        new_dic["absolutePosition"] = \
            struct.unpack("<ddddd", dic["absolutePosition"]["values"][0][4:44])
        _msc = dic["microscopeConditions"]
        new_dic["columnConditions"], _ = \
            self._get_dict(_msc['columnConditions']['values'][0][4:])
        new_dic["stageConditions"], _ = \
            self._get_dict(_msc['stageConditions']['values'][0][4:])
        if "ScanCalibration" in dic:
            new_dic["ScanCalibration"] = \
                struct.unpack("<dd", dic["ScanCalibration"]['values'][0][4:])

        if "region.mapArea" in dic:
            _v0 = struct.unpack("<ldddddd",dic["region.mapArea"]['values'][0])
            # _v1 = struct.unpack("<lldddllllldlll",dic["region.mapArea"]['values'][1])
            new_dic["mapArea"] = _v0
        # print(self._dump_bytes(dic["OINA.Mustang.Data.Regions.region"]["value"][1]))
        # dic["OINA.Mustang.Data.Regions.region"] = \
        #     struct.unpack("<dddd", dic["OINA.Mustang.Data.Regions.region"]["value"][1][4:])
        # self._post_proc_cleanup(dic, ['AreaOfInterest'])
        # new_dic[_uuid] = copy.deepcopy(dic)
        return { _uuid: new_dic }

    def _post_proc_ElectronImage(self, dic):
        print("Dic=",dic)
        _uuid = self._decode_uuid(dic["values"][0][4:20])
        _dat = dic["values"][1]
        if _dat[4] == 0:
            dic['title'], _dat = self._get_string(_dat[5:])
        _s = struct.unpack('<llll', _dat[0:16])
        dic['width'] = _s[1]
        dic['height'] = _s[2]
        dic['offset'] = _s[3]
        dic['_unknown0'] = struct.unpack('<fffff', _dat[16:36])
        dic['File_UUID'] = self._decode_uuid(_dat[36:52])
        dic['_unknown1'] = struct.unpack('<ddd', _dat[56:])
        _dat = dic["imageDataCore.TileSizes.TileSizeInfo"]["values"][0]
        _s = struct.unpack('<llll', _dat[0:16])
        dic['preview_width'] = _s[1]
        dic['preview_height'] = _s[2]
        dic['preview_offset'] = _s[3]
        dic['_unknown_preview'] = self._dump_bytes_debug(_dat[16:])
        _dat = dic["values"][2]
        dic['_unknown2'] = struct.unpack('<llfll', _dat[0:20])
        dic['signal'], _ = self._get_string(_dat[21:])
        dic['pixel'] = dic["ElectronImage.imageDataCore.type"]
        self._post_proc_cleanup(dic, ['values',])
        return { _uuid: dic }



    def __init__(self):
        self.debug = 2

    def decode_metadata(self, cur, table):
        """
        Parameters:
        ----------------
          cur: sqlite3 cursor
          columns: dict
             column type list in record

        Returns:
        ----------------
          metadata: dict
        """

        # SQLite3 DB record structures
        # Data Name, [Supported (1: Fullly, 0: Not, -1: Partly), Column names in DB]
        _DB_TABLES = {
            "Project": [ -1, 'UUID', 'Unknown', 'ProjectTitle', 'Unknown2', 'Data' ],
            "DataNodes": [ 1, 'UUID', 'Unknown', 'Sort', 'UUID_parent',
                           'UUID_data', 'Type', 'Data' ],
            "ExtraData": [ 0, 'UUID', 'Unknown', 'Data' ],
            "Samples": [ 0, 'UUID', 'Unknown', 'SampleName', 'Data' ],
            "Elements": [ 0, 'UUID', 'Unknown', 'Data' ],
            "Analyses": [ 0, 'UUID', 'Unknown', 'AreaName', 'Data' ],
            "ElectronImages": [ 1, 'UUID', 'Unknown', 'ImageName', 'DateTime', 'Data' ],
            "Display256LinearARGBData": [ 0, 'UUID', 'Unknown', 'Data' ],
            "AreaOfInterests": [ 1, 'UUID', 'Unknown', 'Data' ],
            "FsdAcquires": [ 0, 'UUID', 'Unknown', 'DataName', 'Data' ],
            "FsdChannels": [ 0, 'UUID', 'Unknown', 'ChannelName', 'Data' ],
            "FsdMixedImages": [ 0, 'UUID', 'Unknown', 'ImageName', 'Data' ],
            "DisplayARGBData": [0, 'UUID', 'Unknown', 'Data' ],
            "EDXSpectra": [ -1, 'UUID', 'Unknown', 'SpectrumName', 'DateTime', 'Data' ],
            "IdentifiedElements": [ 0, 'UUID', 'Unknown', 'Data' ],
            "QuantResults": [ 0, 'UUID', 'Unknown', 'Data' ],
            "XrayMapImages": [ 0, 'UUID', 'Unknown', 'Xray' ,'DateTime', 'Data' ],
            "LayerMaps": [ 0, 'UUID', 'Unknown', 'Data' ],
            "TileGroups": [ 0, 'UUID', 'Unknown', 'UUID-Unknwon', 'Data' ],
            "FileGroups": [ 0, 'UUID', 'Unknown', '_Data' ],
            "SmartMaps": [ 0, 'UUID', 'Unknown', 'DataName', 'Data', 'Unknown2' ],
            "EbsdMapImages": [ 0, 'UUID', 'Unknown', 'ImageName', 'DateTime', 'Data' ],
            "DisplayIndexedARGBData": [ 0, 'UUID', 'Unknown', 'Data' ],
            "EbsdMaps": [0, 'UUID', 'Unknown', 'DataName', 'Data'],
            "MapAcquires": [0, 'UUID', 'Unknown', 'DataName', 'Data'],
        }
        data = tools.DTBox(box_dots=True, default_box=True)

        # Get data column number from table definition
        if table not in _DB_TABLES:
            if self.debug > 1:
                _logger.warning(f"{tbl} is not found in known DB tables.")
                return None
        column_table = {}
        _columns = _DB_TABLES[table].copy()
        supported = _columns.pop(0)
        for _i, _column in enumerate(_columns):
            column_table[_column] = _i

        if supported == 0:
            if self.debug > 1:
                _logger.warning(f"{table} is not supported yet.")
            return None
        cur.execute("SELECT * FROM " + table)
        _items = []
        _metadata = []
        for row in cur:
            _metadata.append(self.decode_metadata_record(table, row, column_table))
            # remove binary data
            row = list(row)
            del row[column_table["Data"]]
            _items.append(row)

        data[table]["columns"] = column_table.copy()
        data[table]["db_items"] = copy.deepcopy(_items)
        data[table]["metadata"] = copy.deepcopy(_metadata)
        return data

    def decode_metadata_record(self, table, row, column_list):
        bin_dec = _BinaryMetadata()
        _dict = bin_dec.decode_metadata(row, column_list)
        _dict = tools.DTBox(self._convert_metadata_to_dict(_dict),
            box_dots=True, default_box=True)
        if table in self.POST_PROC:
            metadata = self.POST_PROC[table](self, _dict)
        else:
            # currently unsupported tables
            metadata = _dict
            # metadata = {}
        return metadata

    # list of table specific decoders
    POST_PROC = {
        "Project": _post_proc_Project,
        "DataNodes": _post_proc_DataNode,
        "AreaOfInterests": _post_proc_AreaOfInterest,
        "ElectronImages": _post_proc_ElectronImage,
        "XrayMapImages": _post_proc_ElectronImage,
    }



class OxfordOipx:
    def _list_data_nodes(self, data_nodes, _node_uuid, level = 0):
        _node = data_nodes[_node_uuid]
        print("       "[0:level+2], _node["Type"], _node_uuid, end="")
        if "Text" in _node:
            print(_node["Text"], end=", ")
            if "Text2" in _node:
                print(_node["Text2"], end=", ")
            print(" ", _node["UUID_data"])
        else:
            print("<blank>")
        if "children" in _node:
            for _uuid in _node["children"]:
                self._list_data_nodes(data_nodes, data_nodes[_uuid]["UUID"], level + 2)

    def _search_data_nodes(self, data_nodes, _node_uuid):
        for _k, _v in data_nodes.items():
            if _v["UUID_data"] == _node_uuid:
                return data_nodes[_v["UUID_parent"]]["children"]
        return None
    def _get_image_signal_type(self, data_nodes, _node_uuid):
        _children = self._search_data_nodes(data_nodes, _node_uuid)
        if _children is None:
            return None
        for _child in _children:
            if data_nodes[_child]["Type"] == 3:
                return data_nodes[_child]["Text2"]
        return None

    def _decode_data_node(self, cur):
        _data_nodes = []
        top_node = None
        cur.execute("SELECT * FROM DataNodes")
        for i, row in enumerate(cur):
            _uuid, _id1, _sort, _uuid_parent, _uuid_next, _id2, _dat = row
            _dic = self._decode_metadata(_dat)
            _dic = self._post_proc_DataNode(_dic)
            _dic["id1"] = _id1
            _dic["id2"] = _id2
            _dic["sort"] = _sort
            _dic["children"] = []
            _data_nodes.append(copy.deepcopy(_dic))
        _data_nodes = sorted(_data_nodes, key = lambda x: x["sort"])

        data_nodes = {}
        for _nd in _data_nodes:
            data_nodes[_nd['UUID']] = _nd
            _p = _nd['UUID_parent']
            if _p != "00000000-0000-0000-0000-000000000000":
                data_nodes[_p]['children'].append(_nd['UUID'])

        top_node_uuid = _data_nodes[0]['UUID']

        self._list_data_nodes(data_nodes, top_node_uuid, 0)
        return data_nodes, top_node_uuid


    def read_oipx(self, fname):
        data_dir = Path(fname).parent / "data"

        images = []
        data = tools.DTBox(box_dots=True, default_box=True)
        metadata = tools.DTBox()
        metadata_decoder = MetaData()

        # connect to SQLite3 DB
        datanodes = None
        conn = sqlite3.connect(fname)
        if conn is None:
            return None
        cur = conn.cursor()

        # get table list
        cur.execute("SELECT * FROM sqlite_master")
        tables = []
        for row in cur:
            if row[0] == 'table':
                tables.append(row[1])

        _metadata = []
        for tbl in tables:
            _metadata.append(metadata_decoder.decode_metadata(cur, tbl))
        print(_metadata)
        sys.exit(1)
        for tbl in ("ElectronImages", "EDXSpectra"):
            if tbl in self.POST_PROC:
                _proc = (self.POST_PROC)[tbl][1]
                images.extend(_proc(self, metadata[_uuid], row))
        conn.close()
        print("Dataset = ", data)
        return images

    def _set(self, dct, tag_name, val):
        _d = dct
        _x_tags = tag_name.split('.')
        for _tag in _x_tags[:-1]:
            if _tag not in dct:
                dct[_tag] = {}
            dct = dct[_tag]
        _key = _x_tags[-1]
        if _key in dct:
            v0 = dct[_key]
            dct[_key] = [v0]
            dct[_key].append(val)
        else:
            dct[_key] = val
        dct = val

    def _get(self, dct, tag_name):
        _d = dct
        _x_tags = tag_name.split('.')
        for _tag in _x_tags:
            if _tag not in dct:
                return None
            dct = dct[_tag]
        return dct

    def dump_dict(self, dct, lvl=0):
        SPC = "                                              "
        if type(dct) == list:
            for i, itm in enumerate(dct):
                print(SPC[0:lvl*4], f'[{i}] = ')
                self.dump_dict(itm)
        elif type(dct) == dict:
            for _key in dct:
                print(SPC[0:lvl*4], _key, " = ")
                self.dump_dict(dct[_key], lvl + 1)
        else:
            print(SPC[0:lvl*4], dct)


    def _add_block(self, dic, namespace, block_name, block):
        key = ".".join((namespace,block_name)).strip('.')
        dic.add_node(key)
        dic[key] = copy.deepcopy(block)


    def __init__(self):
        pass

    def _encode_uuid(self, text_uuid):
        _t = text_uuid.replace('-', '')
        _uuid = b''
        for i in [3,2,1,0, 5,4, 7,6, 8,9, 10,11,12,13,14,15]:
            _uuid += bytes.fromhex(_t[i*2]+_t[i*2+1])
        return _uuid


    def _read_ElectronImage(self, original_metadata=None, row=None, **kwargs):
        # print(original_metadata)
        images = []
        DATA_TYPE_TABLE = {
            "a:FoldedImageShort": np.uint16,
            "a:FoldedImageFloat": np.float32,
        }

        _om = original_metadata
        _typ = _om['pixel']
        if _typ not in DATA_TYPE_TABLE:
            print("Unknown data type : ", _typ)
            return None
        _dtyp = DATA_TYPE_TABLE[_typ]
        _pix_depth = np.dtype(_dtyp).itemsize
        _uuid = _om["UUID"]
        _w = _om["width"]
        _h = _om["height"]
        _offset = _om["offset"] * _pix_depth  # may be broken when EDS
        _filename = data_dir / (original_metadata["File_UUID"] + ".dat").lower()
        with open(_filename, "rb") as f:
            print(_offset)
            #  f.seek(_offset)  # may be broken when EDS
            data = np.fromfile(f, dtype=_dtyp, count=_w*_h).reshape((_h,_w))
        time_zone = timezone(timedelta(minutes = self.time_zone))
        tm = datetime.strptime(row[3],"%m/%d/%Y %H:%M:%S").replace(tzinfo=time_zone)
        original_metadata["time"] = tm.isoformat()
        _time = tm.time().isoformat()
        _date = tm.date().isoformat()
        _column = _om["AreaOfInterest.columnConditions"]
        _stage = _om["AreaOfInterest.stageConditions"]
        _acq_mode = "SEM" if _column["HighVoltage"] <= 30 else "TEM"
        _om["SigType"] = self._get_image_signal_type(self.datanodes, _uuid)
        _cal = _om["AreaOfInterest.ScanCalibration"][1]
        _mag = _column["Magnification"]
        _fov = _cal / _mag
        _scale = _fov/_w
        _units = ["fm", "pm", "nm", "um", "mm"]
        _unit = _units.pop()
        while _scale < 0.01 and len(_units) > 0:
            _scale *= 1000
            _unit = _units.pop()
        _axes = [
            {
                "name": "y",
                "size": _h,
                "offset": 0.,
                "scale": _scale,
                "units": _unit,
            }, {
                "name": "x",
                "size": _w,
                "offset": 0.,
                "scale": _scale,
                "units": _unit,
            }]
        _metadata = {
            "Acquisition_instrument": {
                _acq_mode: {
                    "Stage": {
                        # "tilt_alpha" : 0,
                        "tilt_beta" : _stage["StageTilt"],
                        "x": _stage["StageX"],
                        "y": _stage["StageY"],
                        "z": _stage["StageZ"],
                    }
                },
                "acquisition_mode": _acq_mode,
                # "beam_current": 0,
                "beam_energy": _column["HighVoltage"],
                "magnification": _mag,
                "microscope": "",
                "working_distance": _column['WorkingDistance'],
            },
            "General": {
                "title": original_metadata["title"],
                "date": _date,
                "time": _time,
                "original_filename": _filename.name,
            },
            "Signal": {
            },
        }
        _dictionary = {
            "data": data,
            "axes": _axes,
            "metadata": _metadata,
            "original_metadata": original_metadata
        }
        images.append(copy.deepcopy(_dictionary))
        return images



    def _read_ebsp(self, filename, **kwargs):
        print("File=", filename)
        with open(filename, "rb") as f:
            _magic = f.read(9)
            if _magic != b"\xfc\xff\xff\xff\xff\xff\xff\xff\x00":
                return None
            toc = []
            _fpos = 9
            _data_pos = os.path.getsize(filename)  # max length
            while True:
                pos = np.fromfile(f, dtype=np.uint64, count=1)
                _fpos += 8
                if pos.size == 0: # broken file
                    break
                pos = pos[0].item() # convert from np.uint64 to int
                if pos < _data_pos:
                    _data_pos = pos
                toc.append(pos)
                if _fpos >= _data_pos:
                    break
            pos = toc[0]
            f.seek(pos)
            _, height, width, _ = np.fromfile(f, dtype=np.uint32, count=4)
            depth = len(toc)
            data = np.ndarray(width * height * depth, dtype=np.uint8)

            for z, pos in enumerate(toc):
                f.seek(16+pos)
                _layer = np.fromfile(f, dtype=np.uint8, count=width * height)
                data[z * width * height:(z + 1) * width * height] = _layer

        _axes = [
            {
                "name": "y",
                "size": height,
                "offset": 0.,
                "scale": 1.,
                "units": "pixels",
            }, {
                "name": "x",
                "size": width,
                "offset": 0.,
                "scale": 1.,
                "units": "pixels",
            }, {
                "name": "z",
                "size": depth,
                "offset": 0.,
                "scale": 1.,
                "units": "pixels",
            }
        ]
        _dictionary = {
            "data": data,
            "axes": _axes,
            "metadata": {},
            "original_metadata": {}
        }
        return _dictionary





def _gen_files(dirname, uuid, exts):
    _files = []
    for ext in exts:
        _file = Path(dirname) / (uuid.lower() + ext)
        if _file.exists():
            _files.append(_file)
    return _files

def file_reader(filename, **kwargs):
    """
    File reader for JEOL Analysist Station software format.

    Parameters
    ----------
    %s
    %s
    rebin_energy : int, Default=1.
        Factor used to rebin the energy dimension. It must be a
        divisor of the number of channels, typically 4096.
    sum_frames : bool, Default=True.
        If ``False``, each individual frame (sweep in JEOL software jargon)
        is loaded. Be aware that loading each individual frame will use a lot of memory.
        However, it can be used in combination with ``rebin_energy``, ``cutoff_at_kV``
        and ``downsample`` to reduce memory usage.
    SI_dtype : dtype, Default=np.uint8
        Set ``dtype`` of the eds dataset. Useful to adjust memory usage
        and maximum number of X-rays per channel.
    cutoff_at_kV : int, float, or None, default=None
        If set (>= 0), use to crop the energy range up to the specified energy.
        If ``None``, the whole energy range is loaded.
        Useful to reduce memory usage.
    downsample : int, Default=1
g        The downsample ratio of the navigation dimension of an EDS
        dataset. It can be an integer or a tuple of length 2 to define ``x`` and ``y``
        separetely and it must be a divisor of the size of the navigation dimension.
    only_valid_data : bool, Default=True
        For ``.pts`` files only. Ignore incomplete and partly
        acquired last frame, which typically occurs when the acquisition was
        interrupted. When loading incomplete data (``only_valid_data=False``),
        the missing data are filled with zeros. If ``sum_frames=True``, this argument
        will be ignored to enforce consistent summing over the mapped area.
    read_em_image : bool, Default=False
        For ``.pts`` files only. If ``True``,
        read SEM/STEM image from ``.pts`` file if available. In this case, both
        the spectrum Image and SEM/STEM Image will be returned as list.
    frame_list : list of integer or None, Default=None
        For ``.pts`` files only. Frames in ``frame_list`` will be loaded.
        For example, ``frame_list=[1,3]`` means second and forth frame will be loaded.
        If ``None``, all frames are loaded.
    frame_shifts : list of [int, int], list of [int, int, int], or None, Default=None
        For ``.pts`` files only. Each frame will be loaded with offset of
        [dy, dx (, and optionary dEnergy)]. Units are pixels/channels.
        The result of estimate_shift2D() can be used as a parameter of frame_shifts.
        This is useful for express drift correction. Not suitable for accurate analysis.

    %s
    """
    file_ext = os.path.splitext(filename)[-1][1:].lower()
    if file_ext not in extension_to_reader_mapping:
        _logger.info(f"{filename} : File type {file_ext} is not supported. Skipping")
        return []
    return extension_to_reader_mapping[file_ext](filename, **kwargs)

file_reader.__doc__ %= (FILENAME_DOC, LAZY_DOC, RETURNS_DOC)

'''
    axes = [
        {
            "name": "Energy",
            "size": header["NumCH"],
            "offset": header["CoefB"],
            "scale": header["CoefA"],
            "units": "keV",
        }
    ]
    metadata = {
        "Acquisition_instrument": {
            mode: {
                "beam_energy": hv,
                "Detector": {
                    "EDS": {
                        "azimuth_angle": footer["Parameters"]["DirAng"],
                        "detector_type": footer["Parameters"]["DetT"],
                        "elevation_angle": footer["Parameters"]["ElevAng"],
                        "energy_resolution_MnKa": footer["Parameters"]["MnKaRES"],
                        "live_time": header["live time"],
                    },
                },
            },
        },
        "General": {
            "original_filename": os.path.basename(filename),
            "date": header["filedate"].date().isoformat(),
            "time": header["filedate"].time().isoformat(),
            "title": "EDX",
        },
        "Signal": {
            "record_by": "spectrum",
            "quantity": "X-rays (Counts)",
            "signal_type": "EDS_" + mode,
        },
    }

    dictionary = {
        "data": data,
        "axes": axes,
        "metadata": metadata,
        "original_metadata": {"Header": header, "Footer": footer},
    }

    return [dictionary]
'''

def _read_oipx(filename, **kwargs):
    oipx = OxfordOipx()
    print("fname=",filename)
    return oipx.read_oipx(filename)

extension_to_reader_mapping = {
    "oipx": _read_oipx,
}
