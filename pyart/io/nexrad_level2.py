"""
Functions for reading NEXRAD level 2 files.
"""

# This file is part of the Py-ART, the Python ARM Radar Toolkit
# https://github.com/ARM-DOE/pyart

# Care has been taken to keep this file free from extraneous dependancies
# so that it can be used by other projects with no/minimal modification.

# Please feel free to use this file in other project provided the license
# below is followed. Keeping the above comment lines would also be helpful
# to direct other back to the Py-ART project and the source of this file.


LICENSE = """
Copyright (c) 2013, UChicago Argonne, LLC
All rights reserved.

Copyright 2013 UChicago Argonne, LLC. This software was produced under U.S.
Government contract DE-AC02-06CH11357 for Argonne National Laboratory (ANL),
which is operated by UChicago Argonne, LLC for the U.S. Department of Energy.
The U.S. Government has rights to use, reproduce, and distribute this
software. NEITHER THE GOVERNMENT NOR UCHICAGO ARGONNE, LLC MAKES ANY
WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
SOFTWARE. If software is modified to produce derivative works, such modified
software should be clearly marked, so as not to confuse it with the version
available from ANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions
are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name of UChicago Argonne, LLC, Argonne National
      Laboratory, ANL, the U.S. Government, nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY UCHICAGO ARGONNE, LLC AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL UCHICAGO ARGONNE, LLC OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import bz2
import struct
import warnings
from datetime import datetime, timedelta

import numpy as np

# Waveform type codes, MSG_5 per-cut data (E2 lower byte).
# ICD 2620002AA (Build 24.0), Table XI and Table XVIII, footnote (1) on
# Table XVIII: type 3 uses the PRFs of type 2, type 4 (Batch) uses the
# type 1 PRFs for its surveillance part and the type 2 PRFs for its Doppler part.
_WAVEFORM_TYPES = {
    0: "not_applicable",
    1: "contiguous_surveillance",  # CS
    2: "contiguous_doppler_with_ambiguity_resolution",  # CD/W
    3: "contiguous_doppler_without_ambiguity_resolution",  # CD/WO
    4: "batch",  # B
    5: "staggered_pulse_pair",
}


# Decoders for the SAILS/MRLE/MPDA/base-tilt bit fields.
# ICD 2620002AA, Table XI and notes (14)-(17). Bit 0 is the least significant
# bit (note (4) of Table XI).
def _decode_vcp_sequencing(value):
    """
    Decode the VCP sequencing halfword (MSG_5 HW 9, Table XI note 15).

    This field is reserved for RPG use (note 14). It is zero unless the VCP is
    part of an active VCP sequence, so it is zero in most files. Do not rely on
    it for the number of cuts in a volume, use MSG_5 HW 4 (num_cuts) instead.
    """
    return {
        "num_elevations": value & 0x1F,  # bits 0-4 (truncated VCPs)
        "max_sails_cuts": (value >> 5) & 0x03,  # bits 5-6
        "sequence_active": bool(value & 0x2000),  # bit 13
        "truncated_vcp": bool(value & 0x4000),  # bit 14
    }


def _decode_vcp_supplemental(value):
    """
    Decode the VCP supplemental halfword (MSG_5 HW 10, Table XI note 16).

    Bits 8-10 are spares according to note 16 (the table labels bit 10 as
    "MPDA Cuts Added"), so they are deliberately not decoded.
    SAILS and MRLE cannot be active at the same time.
    """
    return {
        "sails_vcp": bool(value & 0x0001),  # bit 0
        "num_sails_cuts": (value >> 1) & 0x07,  # bits 1-3 (max 3)
        "mrle_vcp": bool(value & 0x0010),  # bit 4
        "num_mrle_cuts": (value >> 5) & 0x07,  # bits 5-7 (max 4)
        "mpda_vcp": bool(value & 0x0800),  # bit 11
        "base_tilt_vcp": bool(value & 0x1000),  # bit 12
        "num_base_tilts": (value >> 13) & 0x07,  # bits 13-15
    }


def _decode_elevation_supplemental(value):
    """
    Decode the per-cut supplemental halfword (MSG_5 E15, Table XI note 17).

    A cut cannot be both a SAILS and an MRLE cut. For a SAILS cut, bits 1-3 are
    the SAILS sequence number; split (CS/CD) cuts share one number. For an MRLE
    cut, bits 5-7 equal the RPG elevation index of the cut, which is not a
    1..N counter.
    """
    return {
        "sails_cut": bool(value & 0x0001),  # bit 0
        "sails_sequence_number": (value >> 1) & 0x07,  # bits 1-3
        "mrle_cut": bool(value & 0x0010),  # bit 4
        "mrle_sequence_number": (value >> 5) & 0x07,  # bits 5-7
        "mpda_cut": bool(value & 0x0200),  # bit 9 (bit 8 is a spare)
        "base_tilt_cut": bool(value & 0x0400),  # bit 10
    }


class NEXRADLevel2File:
    """
    Class for accessing data in a NEXRAD (WSR-88D) Level II file.

    NEXRAD Level II files [1]_, also know as NEXRAD Archive Level II or
    WSR-88D Archive level 2, are available from the NOAA National Climate Data
    Center [2]_ as well as on the UCAR THREDDS Data Server [3]_. Files with
    uncompressed messages and compressed messages are supported. This class
    supports reading both "message 31" and "message 1" type files.

    Parameters
    ----------
    filename : str
        Filename of Archive II file to read.

    Attributes
    ----------
    radial_records : list
        Radial (1 or 31) messages in the file.
    nscans : int
        Number of scans in the file.
    scan_msgs : list of arrays
        Each element specifies the indices of the message in the
        radial_records attribute which belong to a given scan.
    volume_header : dict
        Volume header.
    vcp : dict
        VCP information dictionary.
    _records : list
        A list of all records (message) in the file.
    _fh : file-like
        File like object from which data is read.
    _msg_type : '31' or '1':
        Type of radial messages in file.

    References
    ----------
    .. [1] http://www.roc.noaa.gov/WSR88D/Level_II/Level2Info.aspx
    .. [2] http://www.ncdc.noaa.gov/
    .. [3] http://thredds.ucar.edu/thredds/catalog.html

    """

    def __init__(self, filename):
        """initalize the object."""
        # read in the volume header and compression_record
        if hasattr(filename, "read"):
            fh = filename
        else:
            fh = open(filename, "rb")
        size = _structure_size(VOLUME_HEADER)
        self.volume_header = _unpack_structure(fh.read(size), VOLUME_HEADER)
        compression_record = fh.read(COMPRESSION_RECORD_SIZE)

        # read the records in the file, decompressing as needed
        compression_slice = slice(CONTROL_WORD_SIZE, CONTROL_WORD_SIZE + 2)
        compression_or_ctm_info = compression_record[compression_slice]
        if compression_or_ctm_info == b"BZ":
            buf = _decompress_records(fh)
        # The 12-byte compression record previously held the Channel Terminal
        # Manager (CTM) information. Bytes 4 through 6 contain the size of the
        # record (2432) as a big endian unsigned short, which is encoded as
        # b'\t\x80' == struct.pack('>H', 2432).
        # Newer files zero out this section.
        elif compression_or_ctm_info in (b"\x00\x00", b"\t\x80"):
            buf = fh.read()
        else:
            raise OSError("unknown compression record")
        self._fh = fh

        # read the records from the buffer
        self._records = []
        buf_length = len(buf)
        pos = 0
        while pos < buf_length:
            pos, dic = _get_record_from_buf(buf, pos)
            self._records.append(dic)

        # pull out radial records (1 or 31) which contain the moment data.
        self.radial_records = [r for r in self._records if r["header"]["type"] == 31]
        self._msg_type = "31"
        if len(self.radial_records) == 0:
            self.radial_records = [r for r in self._records if r["header"]["type"] == 1]
            self._msg_type = "1"
        if len(self.radial_records) == 0:
            raise ValueError("No MSG31 records found, cannot read file")
        elev_nums = np.array(
            [m["msg_header"]["elevation_number"] for m in self.radial_records]
        )
        self.scan_msgs = [
            np.where(elev_nums == i + 1)[0] for i in range(elev_nums.max())
        ]
        self.nscans = len(self.scan_msgs)

        # pull out the vcp record
        msg_5 = [r for r in self._records if r["header"]["type"] == 5]

        if len(msg_5):
            self.vcp = msg_5[0]
        else:
            # There is no VCP Data.. This is uber dodgy
            warnings.warn(
                "No MSG5 detected. Setting to meaningless data. "
                "Rethink your life choices and be ready for errors."
                "Specifically fixed angle data will be missing"
            )

            self.vcp = None
        return

    def close(self):
        """Close the file."""
        self._fh.close()

    def location(self):
        """
        Find the location of the radar.

        Returns all zeros if location is not available.

        Returns
        -------
        latitude : float
            Latitude of the radar in degrees.
        longitude : float
            Longitude of the radar in degrees.
        height : int
            Height of radar and feedhorn in meters above mean sea level.

        """
        if self._msg_type == "31":
            dic = self.radial_records[0]["VOL"]
            height = dic["height"] + dic["feedhorn_height"]
            return dic["lat"], dic["lon"], height
        else:
            return 0.0, 0.0, 0.0

    def scan_info(self, scans=None):
        """
        Return a list of dictionaries with scan information.

        Parameters
        ----------
        scans : list ot None
            Scans (0 based) for which ray (radial) azimuth angles will be
            retrieved.  None (the default) will return the angles for all
            scans in the volume.

        Returns
        -------
        scan_info : list, optional
            A list of the scan performed with a dictionary with keys
            'moments', 'ngates', 'nrays', 'first_gate' and 'gate_spacing'
            for each scan.  The 'moments', 'ngates', 'first_gate', and
            'gate_spacing' keys are lists of the NEXRAD moments and gate
            information for that moment collected during the specific scan.
            The 'nrays' key provides the number of radials collected in the
            given scan.

        """
        info = []

        if scans is None:
            scans = range(self.nscans)
        for scan in scans:
            nrays = self.get_nrays(scan)
            if nrays < 2:
                # NOTE: this mutates self.nscans but not scan_msgs. It is left
                # unchanged because existing readers may depend on it; the new
                # get_scan_types() deliberately does not use self.nscans.
                self.nscans -= 1
                continue
            msg31_number = self.scan_msgs[scan][0]
            msg = self.radial_records[msg31_number]
            nexrad_moments = ["REF", "VEL", "SW", "ZDR", "PHI", "RHO", "CFP"]
            moments = [f for f in nexrad_moments if f in msg]
            ngates = [msg[f]["ngates"] for f in moments]
            gate_spacing = [msg[f]["gate_spacing"] for f in moments]
            first_gate = [msg[f]["first_gate"] for f in moments]
            info.append(
                {
                    "nrays": nrays,
                    "ngates": ngates,
                    "gate_spacing": gate_spacing,
                    "first_gate": first_gate,
                    "moments": moments,
                }
            )
        return info

    def get_vcp_pattern(self):
        """
        Return the numerical volume coverage pattern (VCP) or None if unknown.
        """
        if self.vcp is None:
            return None
        else:
            return self.vcp["msg5_header"]["pattern_number"]

    def get_vcp_flags(self):
        """
        Return decoded VCP-level SAILS/MRLE/MPDA/base-tilt information.

        Returns
        -------
        flags : dict or None
            Dictionary with keys 'sequencing' and 'supplemental', each a dict
            of decoded values (see _decode_vcp_sequencing and
            _decode_vcp_supplemental). None if no MSG_5 header is available.

        """
        if self.vcp is None or "msg5_header" not in self.vcp:
            return None
        hdr = self.vcp["msg5_header"]
        return {
            "sequencing": _decode_vcp_sequencing(hdr["vcp_sequencing"]),
            "supplemental": _decode_vcp_supplemental(hdr["vcp_supplemental"]),
        }

    def get_scan_types(self, scans=None, elevation_tolerance=1.0):
        """
        Return per-scan metadata: waveform type and SAILS/MRLE flags.

        Scans are matched to MSG_5 cuts by the radial's elevation_number (not
        by list position). If the MSG_5 nominal elevation differs from the
        measured elevation by more than ``elevation_tolerance`` degrees, the
        MSG_5 information is not attached for that scan (waveform_type and
        waveform_code stay None) and a RuntimeWarning is issued.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) to describe. None (the default) uses all scans in
            the volume, including scans that scan_info() may have skipped.
        elevation_tolerance : float
            Maximum allowed difference in degrees between the MSG_5 nominal
            elevation and the measured median elevation.

        Returns
        -------
        scan_types : list of dict
            One entry per scan with keys 'scan', 'elevation_number',
            'elevation', 'waveform_type', 'waveform_code', 'sails_cut',
            'sails_sequence_number', 'mrle_cut', 'mrle_sequence_number',
            'mpda_cut' and 'base_tilt_cut'. 'waveform_type' is a name from
            _WAVEFORM_TYPES and 'waveform_code' the raw ICD code (1 = CS,
            2 = CD/W, 3 = CD/WO, 4 = Batch, 5 = Staggered Pulse Pair).
            Scans with no radials (gaps in the elevation numbers) are skipped.

        """
        if scans is None:
            scans = range(len(self.scan_msgs))  # not self.nscans
        cuts = []
        if self.vcp is not None:
            cuts = self.vcp.get("cut_parameters", [])
        angle_scale = 1.0 if self._msg_type == "31" else 180 / (4096 * 8.0)

        out = []
        for scan in scans:
            msg_idx = self.scan_msgs[scan]
            if len(msg_idx) == 0:  # gap in elevation numbers
                continue
            hdr = self.radial_records[msg_idx[0]]["msg_header"]
            elev_num = int(hdr["elevation_number"])
            measured = float(
                np.median(
                    [
                        self.radial_records[i]["msg_header"]["elevation_angle"]
                        for i in msg_idx
                    ]
                )
                * angle_scale
            )

            entry = {
                "scan": scan,
                "elevation_number": elev_num,
                "elevation": measured,
                "waveform_type": None,
                "waveform_code": None,
                "sails_cut": None,
                "sails_sequence_number": None,
                "mrle_cut": None,
                "mrle_sequence_number": None,
                "mpda_cut": None,
                "base_tilt_cut": None,
            }

            # MSG_5 cut data only exists for MSG31 files
            cut_idx = elev_num - 1
            if self._msg_type == "31" and 0 <= cut_idx < len(cuts):
                cut = cuts[cut_idx]
                nominal = cut["elevation_angle"] * 360.0 / 65536.0
                if abs(nominal - measured) <= elevation_tolerance:
                    wf = cut["waveform_type"]
                    entry["waveform_code"] = int(wf)
                    entry["waveform_type"] = _WAVEFORM_TYPES.get(wf, str(wf))
                    supplemental = cut.get("supplemental_data")
                    if supplemental is not None:
                        entry.update(_decode_elevation_supplemental(supplemental))
                else:
                    warnings.warn(
                        f"Scan {scan}: MSG_5 cut {elev_num} nominal elevation "
                        f"{nominal:.2f} deg does not match measured "
                        f"{measured:.2f} deg; cut metadata not attached.",
                        RuntimeWarning,
                    )
            out.append(entry)
        return out

    def get_scans_by_type(self, elevation=None, waveform=None, atol=0.25, **flags):
        """
        Return indices of scans matching elevation, waveform type and/or flags.

        Parameters
        ----------
        elevation : float or None
            Elevation in degrees to match within ``atol``. None matches all.
            Some sites have a supplemental low elevation (e.g. 0.3 or 0.4
            degrees, ICD 2620002AA Appendix C), adjust ``atol`` if needed.
        waveform : str, iterable of str, or None
            One or several values of _WAVEFORM_TYPES, e.g.
            'contiguous_surveillance'. Contiguous Doppler is two types
            (with and without ambiguity resolution), pass both to match
            either. None matches all.
        atol : float
            Elevation tolerance in degrees.
        **flags
            Additional exact matches on get_scan_types() keys, e.g.
            ``sails_cut=False``.

        Returns
        -------
        scans : list of int
            Scan indices (0 based), usable as the ``scans`` argument of
            get_data, get_elevation_angles, etc.

        Examples
        --------
        >>> nf.get_scans_by_type(0.5, "contiguous_surveillance")
        >>> nf.get_scans_by_type(0.5, sails_cut=True)
        >>> nf.get_scans_by_type(
        ...     0.5,
        ...     [
        ...         "contiguous_doppler_with_ambiguity_resolution",
        ...         "contiguous_doppler_without_ambiguity_resolution",
        ...     ],
        ... )

        """
        allowed = None
        if waveform is not None:
            allowed = {waveform} if isinstance(waveform, str) else set(waveform)

        result = []
        for entry in self.get_scan_types():
            if elevation is not None and abs(entry["elevation"] - elevation) > atol:
                continue
            if allowed is not None and entry["waveform_type"] not in allowed:
                continue
            if any(entry.get(k) != v for k, v in flags.items()):
                continue
            result.append(entry["scan"])
        return result

    def get_nrays(self, scan):
        """
        Return the number of rays in a given scan.

        Parameters
        ----------
        scan : int
            Scan of interest (0 based).

        Returns
        -------
        nrays : int
            Number of rays (radials) in the scan.

        """
        return len(self.scan_msgs[scan])

    def get_range(self, scan_num, moment):
        """
        Return an array of gate ranges for a given scan and moment.

        Parameters
        ----------
        scan_num : int
            Scan number (0 based).
        moment : 'REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO', or 'CFP'
            Moment of interest.

        Returns
        -------
        range : ndarray
            Range in meters from the antenna to the center of gate (bin).

        """
        dic = self.radial_records[self.scan_msgs[scan_num][0]][moment]
        ngates = dic["ngates"]
        first_gate = dic["first_gate"]
        gate_spacing = dic["gate_spacing"]
        return np.arange(ngates) * gate_spacing + first_gate

    # helper functions for looping over scans
    def _msg_nums(self, scans):
        """Find the all message number for a list of scans."""
        return np.concatenate([self.scan_msgs[i] for i in scans])

    def _radial_array(self, scans, key):
        """
        Return an array of radial header elements for all rays in scans.
        """
        msg_nums = self._msg_nums(scans)
        temp = [self.radial_records[i]["msg_header"][key] for i in msg_nums]
        return np.array(temp)

    def _radial_sub_array(self, scans, key):
        """
        Return an array of RAD or msg_header elements for all rays in scans.
        """
        msg_nums = self._msg_nums(scans)
        if self._msg_type == "31":
            tmp = [self.radial_records[i]["RAD"][key] for i in msg_nums]
        else:
            tmp = [self.radial_records[i]["msg_header"][key] for i in msg_nums]
        return np.array(tmp)

    def get_times(self, scans=None):
        """
        Retrieve the times at which the rays were collected.

        Parameters
        ----------
        scans : list or None
            Scans (0-based) to retrieve ray (radial) collection times from.
            None (the default) will return the times for all scans in the
            volume.

        Returns
        -------
        time_start : Datetime
            Initial time.
        time : ndarray
            Offset in seconds from the initial time at which the rays
            in the requested scans were collected.

        """
        if scans is None:
            scans = range(self.nscans)
        days = self._radial_array(scans, "collect_date")
        secs = self._radial_array(scans, "collect_ms") / 1000.0
        offset = timedelta(days=int(days[0]) - 1, seconds=int(secs[0]))
        time_start = datetime(1970, 1, 1) + offset
        time = secs - int(secs[0]) + (days - days[0]) * 86400
        return time_start, time

    def get_azimuth_angles(self, scans=None):
        """
        Retrieve the azimuth angles of all rays in the requested scans.

        Parameters
        ----------
        scans : list ot None
            Scans (0 based) for which ray (radial) azimuth angles will be
            retrieved. None (the default) will return the angles for all
            scans in the volume.

        Returns
        -------
        angles : ndarray
            Azimuth angles in degress for all rays in the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        if self._msg_type == "1":
            scale = 180 / (4096 * 8.0)
        else:
            scale = 1.0
        return self._radial_array(scans, "azimuth_angle") * scale

    def get_elevation_angles(self, scans=None):
        """
        Retrieve the elevation angles of all rays in the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which ray (radial) azimuth angles will be
            retrieved. None (the default) will return the angles for
            all scans in the volume.

        Returns
        -------
        angles : ndarray
            Elevation angles in degress for all rays in the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        if self._msg_type == "1":
            scale = 180 / (4096 * 8.0)
        else:
            scale = 1.0
        return self._radial_array(scans, "elevation_angle") * scale

    def get_target_angles(self, scans=None):
        """
        Retrieve the target elevation angle of the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which the target elevation angles will be
            retrieved. None (the default) will return the angles for all
            scans in the volume.

        Returns
        -------
        angles : ndarray
            Target elevation angles in degress for the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        if self._msg_type == "31":
            if self.vcp is not None:
                cut_parameters = self.vcp["cut_parameters"]
            else:
                cut_parameters = [{"elevation_angle": 0.0}] * self.nscans
            scale = 360.0 / 65536.0
            return np.array(
                [cut_parameters[i]["elevation_angle"] * scale for i in scans],
                dtype="float32",
            )
        else:
            scale = 180 / (4096 * 8.0)
            msgs = [self.radial_records[self.scan_msgs[i][0]] for i in scans]
            return np.round(
                np.array(
                    [m["msg_header"]["elevation_angle"] * scale for m in msgs],
                    dtype="float32",
                ),
                1,
            )

    def get_nyquist_vel(self, scans=None):
        """
        Retrieve the Nyquist velocities of the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which the Nyquist velocities will be
            retrieved. None (the default) will return the velocities for all
            scans in the volume.

        Returns
        -------
        velocities : ndarray
            Nyquist velocities (in m/s) for the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        return self._radial_sub_array(scans, "nyquist_vel") * 0.01

    def get_unambigous_range(self, scans=None):
        """
        Retrieve the unambiguous range of the requested scans.

        Parameters
        ----------
        scans : list or None
            Scans (0 based) for which the unambiguous range will be retrieved.
            None (the default) will return the range for all scans in the
            volume.

        Returns
        -------
        unambiguous_range : ndarray
            Unambiguous range (in meters) for the requested scans.

        """
        if scans is None:
            scans = range(self.nscans)
        # unambiguous range is stored in tenths of km, x100 for meters
        return self._radial_sub_array(scans, "unambig_range") * 100.0

    def get_data(self, moment, max_ngates, scans=None, raw_data=False):
        """
        Retrieve moment data for a given set of scans.

        Masked points indicate that the data was not collected, below
        threshold or is range folded.

        Parameters
        ----------
        moment : 'REF', 'VEL', 'SW', 'ZDR', 'PHI', 'RHO', or 'CFP'
            Moment for which to to retrieve data.
        max_ngates : int
            Maximum number of gates (bins) in any ray.
            requested.
        raw_data : bool
            True to return the raw data, False to perform masking as well as
            applying the appropiate scale and offset to the data.  When
            raw_data is True values of 1 in the data likely indicate that
            the gate was not present in the sweep, in some cases in will
            indicate range folded data.
        scans : list or None.
            Scans to retrieve data from (0 based). None (the default) will
            get the data for all scans in the volume.

        Returns
        -------
        data : ndarray

        """
        if scans is None:
            scans = range(self.nscans)

        # determine the number of rays
        msg_nums = self._msg_nums(scans)
        nrays = len(msg_nums)
        # extract the data
        set_datatype = False
        data = np.ones((nrays, max_ngates), ">B")
        for i, msg_num in enumerate(msg_nums):
            msg = self.radial_records[msg_num]
            if moment not in msg.keys():
                continue
            if not set_datatype:
                data = data.astype(">" + _bits_to_code(msg, moment))
                set_datatype = True

            ngates = min(msg[moment]["ngates"], max_ngates, len(msg[moment]["data"]))
            data[i, :ngates] = msg[moment]["data"][:ngates]
        # return raw data if requested
        if raw_data:
            return data

        # mask, scan and offset, assume that the offset and scale
        # are the same in all scans/gates
        for scan in scans:  # find a scan which contains the moment
            msg_num = self.scan_msgs[scan][0]
            msg = self.radial_records[msg_num]
            if moment in msg.keys():
                offset = np.float32(msg[moment]["offset"])
                scale = np.float32(msg[moment]["scale"])
                mask = data <= 1
                scaled_data = (data - offset) / scale
                return np.ma.array(scaled_data, mask=mask)

        # moment is not present in any scan, mask all values
        return np.ma.masked_less_equal(data, 1)


def _bits_to_code(msg, moment):
    """
    Convert number of bits to the proper code for unpacking.
    Based on the code found in MetPy:
    https://github.com/Unidata/MetPy/blob/40d5c12ab341a449c9398508bd41
    d010165f9eeb/src/metpy/io/_tools.py#L313-L321
    """
    if msg["header"]["type"] == 1:
        word_size = msg[moment]["data"].dtype
        if word_size == "uint16":
            return "H"
        elif word_size == "uint8":
            return "B"
        else:
            warnings.warn(('Unsupported bit size: %s. Returning "B"', word_size))
            return "B"

    elif msg["header"]["type"] == 31:
        word_size = msg[moment]["word_size"]
        if word_size == 16:
            return "H"
        elif word_size == 8:
            return "B"
        else:
            warnings.warn(('Unsupported bit size: %s. Returning "B"', word_size))
            return "B"
    else:
        raise TypeError("Unsupported msg type %s", msg["header"]["type"])


def _decompress_records(file_handler):
    """
    Decompressed the records from an BZ2 compressed Archive 2 file.
    """
    file_handler.seek(0)
    cbuf = file_handler.read()  # read all data from the file
    decompressor = bz2.BZ2Decompressor()
    skip = _structure_size(VOLUME_HEADER) + CONTROL_WORD_SIZE
    buf = bytearray(decompressor.decompress(cbuf[skip:]))
    while len(decompressor.unused_data):
        cbuf = decompressor.unused_data
        decompressor = bz2.BZ2Decompressor()
        buf += decompressor.decompress(cbuf[CONTROL_WORD_SIZE:])

    return buf[COMPRESSION_RECORD_SIZE:]


def _get_record_from_buf(buf, pos):
    """Retrieve and unpack a NEXRAD record from a buffer."""
    dic = {"header": _unpack_from_buf(buf, pos, MSG_HEADER)}
    msg_type = dic["header"]["type"]

    if msg_type == 31:
        new_pos = _get_msg31_from_buf(buf, pos, dic)
    elif msg_type == 5:
        # Sometimes we encounter incomplete buffers
        try:
            new_pos = _get_msg5_from_buf(buf, pos, dic)
        except struct.error:
            warnings.warn(
                "Encountered incomplete MSG5. File may be corrupt.", RuntimeWarning
            )
            new_pos = pos + RECORD_SIZE
    elif msg_type == 29:
        new_pos = _get_msg29_from_buf(pos, dic)
        warnings.warn("Message 29 encountered, not parsing.", RuntimeWarning)
    elif msg_type == 1:
        new_pos = _get_msg1_from_buf(buf, pos, dic)
    else:  # not message 31 or 1, no decoding performed
        new_pos = pos + RECORD_SIZE

    return new_pos, dic


def _get_msg29_from_buf(pos, dic):
    msg_size = dic["header"]["size"]
    if msg_size == 65535:
        msg_size = dic["header"]["segments"] << 16 | dic["header"]["seg_num"]
    msg_header_size = _structure_size(MSG_HEADER)
    new_pos = pos + msg_header_size + msg_size
    return new_pos


def _get_msg31_from_buf(buf, pos, dic):
    """Retrieve and unpack a MSG31 record from a buffer."""
    msg_size = dic["header"]["size"] * 2 - 4
    msg_header_size = _structure_size(MSG_HEADER)
    new_pos = pos + msg_header_size + msg_size
    mbuf = buf[pos + msg_header_size : new_pos]
    msg_31_header = _unpack_from_buf(mbuf, 0, MSG_31)
    block_pointers = [
        v for k, v in msg_31_header.items() if k.startswith("block_pointer") and v > 0
    ]
    for block_pointer in block_pointers:
        block_name, block_dic = _get_msg31_data_block(mbuf, block_pointer)
        dic[block_name] = block_dic

    dic["msg_header"] = msg_31_header
    return new_pos


def _get_msg31_data_block(buf, ptr):
    """Unpack a msg_31 data block into a dictionary."""
    block_name = buf[ptr + 1 : ptr + 4].decode("ascii").strip()

    if block_name == "VOL":
        dic = _unpack_from_buf(buf, ptr, VOLUME_DATA_BLOCK)
    elif block_name == "ELV":
        dic = _unpack_from_buf(buf, ptr, ELEVATION_DATA_BLOCK)
    elif block_name == "RAD":
        dic = _unpack_from_buf(buf, ptr, RADIAL_DATA_BLOCK)
    elif block_name in ["REF", "VEL", "SW", "ZDR", "PHI", "RHO", "CFP"]:
        dic = _unpack_from_buf(buf, ptr, GENERIC_DATA_BLOCK)
        ngates = dic["ngates"]
        ptr2 = ptr + _structure_size(GENERIC_DATA_BLOCK)
        if dic["word_size"] == 16:
            data = np.frombuffer(buf[ptr2 : ptr2 + ngates * 2], ">u2")
        elif dic["word_size"] == 8:
            data = np.frombuffer(buf[ptr2 : ptr2 + ngates], ">u1")
        else:
            warnings.warn(
                'Unsupported bit size: %s. Returning array dtype "B"', dic["word_size"]
            )
        dic["data"] = data
    else:
        dic = {}
    return block_name, dic


def _get_msg1_from_buf(buf, pos, dic):
    """Retrieve and unpack a MSG1 record from a buffer."""
    msg_header_size = _structure_size(MSG_HEADER)
    msg1_header = _unpack_from_buf(buf, pos + msg_header_size, MSG_1)
    dic["msg_header"] = msg1_header

    sur_nbins = int(msg1_header["sur_nbins"])
    doppler_nbins = int(msg1_header["doppler_nbins"])

    sur_step = int(msg1_header["sur_range_step"])
    doppler_step = int(msg1_header["doppler_range_step"])

    sur_first = int(msg1_header["sur_range_first"])
    doppler_first = int(msg1_header["doppler_range_first"])
    if doppler_first > 2**15:
        doppler_first = doppler_first - 2**16

    if msg1_header["sur_pointer"]:
        offset = pos + msg_header_size + msg1_header["sur_pointer"]
        data = np.frombuffer(buf[offset : offset + sur_nbins], ">u1")
        dic["REF"] = {
            "ngates": sur_nbins,
            "gate_spacing": sur_step,
            "first_gate": sur_first,
            "data": data,
            "scale": 2.0,
            "offset": 66.0,
        }
    if msg1_header["vel_pointer"]:
        offset = pos + msg_header_size + msg1_header["vel_pointer"]
        data = np.frombuffer(buf[offset : offset + doppler_nbins], ">u1")
        dic["VEL"] = {
            "ngates": doppler_nbins,
            "gate_spacing": doppler_step,
            "first_gate": doppler_first,
            "data": data,
            "scale": 2.0,
            "offset": 129.0,
        }
        if msg1_header["doppler_resolution"] == 4:
            # 1 m/s resolution velocity, offset remains 129.
            dic["VEL"]["scale"] = 1.0
    if msg1_header["width_pointer"]:
        offset = pos + msg_header_size + msg1_header["width_pointer"]
        data = np.frombuffer(buf[offset : offset + doppler_nbins], ">u1")
        dic["SW"] = {
            "ngates": doppler_nbins,
            "gate_spacing": doppler_step,
            "first_gate": doppler_first,
            "data": data,
            "scale": 2.0,
            "offset": 129.0,
        }
    return pos + RECORD_SIZE


def _get_msg5_from_buf(buf, pos, dic):
    """Retrieve and unpack a MSG5 record from a buffer."""
    msg_header_size = _structure_size(MSG_HEADER)
    msg5_header_size = _structure_size(MSG_5)
    msg5_elev_size = _structure_size(MSG_5_ELEV)

    dic["msg5_header"] = _unpack_from_buf(buf, pos + msg_header_size, MSG_5)
    dic["cut_parameters"] = []
    for i in range(dic["msg5_header"]["num_cuts"]):
        pos2 = pos + msg_header_size + msg5_header_size + msg5_elev_size * i
        dic["cut_parameters"].append(_unpack_from_buf(buf, pos2, MSG_5_ELEV))
    return pos + RECORD_SIZE


def _structure_size(structure):
    """Find the size of a structure in bytes."""
    return struct.calcsize(">" + "".join([i[1] for i in structure]))


def _unpack_from_buf(buf, pos, structure):
    """Unpack a structure from a buffer."""
    size = _structure_size(structure)
    return _unpack_structure(buf[pos : pos + size], structure)


def _unpack_structure(string, structure):
    """Unpack a structure from a string."""
    fmt = ">" + "".join([i[1] for i in structure])  # NEXRAD is big-endian
    lst = struct.unpack(fmt, string)
    return dict(zip([i[0] for i in structure], lst))


# NEXRAD Level II file structures and sizes
# The details on these structures are documented in:
# "Interface Control Document for the Archive II/User", RDA Build 23.0
# Document Number 2620010J
# and
# "Interface Control Document for the RDA/RPG", Build 24.0
# Document Number 2620002AA
# Tables and page number refer to those in the second document unless
# otherwise noted.
RECORD_SIZE = 2432
COMPRESSION_RECORD_SIZE = 12
CONTROL_WORD_SIZE = 4

# format of structure elements
# section 3.2.1, page 3-2
CODE1 = "B"
CODE2 = "H"
INT1 = "B"
INT2 = "H"
INT4 = "I"
REAL4 = "f"
REAL8 = "d"
SINT1 = "b"
SINT2 = "h"
SINT4 = "i"

# Figure 1 in Interface Control Document for the Archive II/User
# page 7-2
VOLUME_HEADER = (
    ("tape", "9s"),
    ("extension", "3s"),
    ("date", "I"),
    ("time", "I"),
    ("icao", "4s"),
)

# Table II Message Header Data
# page 3-7
MSG_HEADER = (
    ("size", INT2),  # size of data, no including header
    ("channels", INT1),
    ("type", INT1),
    ("seq_id", INT2),
    ("date", INT2),
    ("ms", INT4),
    ("segments", INT2),
    ("seg_num", INT2),
)

# Table XVII Digital Radar Generic Format Blocks (Message Type 31)
# Table XVII-A, pages 3-81 to 3-84
MSG_31 = (
    ("id", "4s"),  # 0-3
    ("collect_ms", INT4),  # 4-7
    ("collect_date", INT2),  # 8-9
    ("azimuth_number", INT2),  # 10-11
    ("azimuth_angle", REAL4),  # 12-15
    ("compress_flag", CODE1),  # 16
    ("spare_0", INT1),  # 17
    ("radial_length", INT2),  # 18-19
    ("azimuth_resolution", CODE1),  # 20  1 = 0.5 deg, 2 = 1.0 deg
    # 21  Table III-C: 0 new elevation, 1 intermediate, 2 end of elevation,
    #     3 start of volume, 4 end of volume, 5 new elevation (last in VCP);
    #     add 0x80 for bad data. (Formerly named "radial_spacing".)
    ("radial_status", CODE1),
    ("elevation_number", INT1),  # 22  1 to 32, elevation number within volume
    ("cut_sector", INT1),  # 23
    ("elevation_angle", REAL4),  # 24-27
    ("radial_blanking", CODE1),  # 28
    ("azimuth_mode", SINT1),  # 29
    ("block_count", INT2),  # 30-31
    ("block_pointer_1", INT4),  # 32-35  Volume Data Constant XVII-E
    ("block_pointer_2", INT4),  # 36-39  Elevation Data Constant XVII-F
    ("block_pointer_3", INT4),  # 40-43  Radial Data Constant XVII-H
    ("block_pointer_4", INT4),  # 44-47  Moment "REF" XVII-{B/I}
    ("block_pointer_5", INT4),  # 48-51  Moment "VEL"
    ("block_pointer_6", INT4),  # 52-55  Moment "SW"
    ("block_pointer_7", INT4),  # 56-59  Moment "ZDR"
    ("block_pointer_8", INT4),  # 60-63  Moment "PHI"
    ("block_pointer_9", INT4),  # 64-67  Moment "RHO"
    ("block_pointer_10", INT4),  # 68-71  Moment "CFP"
)


# Table III Digital Radar Data (Message Type 1)
# pages 3-7 to
MSG_1 = (
    ("collect_ms", INT4),  # 0-3
    ("collect_date", INT2),  # 4-5
    ("unambig_range", SINT2),  # 6-7
    ("azimuth_angle", CODE2),  # 8-9
    ("azimuth_number", INT2),  # 10-11
    ("radial_status", CODE2),  # 12-13
    ("elevation_angle", INT2),  # 14-15
    ("elevation_number", INT2),  # 16-17
    ("sur_range_first", CODE2),  # 18-19
    ("doppler_range_first", CODE2),  # 20-21
    ("sur_range_step", CODE2),  # 22-23
    ("doppler_range_step", CODE2),  # 24-25
    ("sur_nbins", INT2),  # 26-27
    ("doppler_nbins", INT2),  # 28-29
    ("cut_sector_num", INT2),  # 30-31
    ("calib_const", REAL4),  # 32-35
    ("sur_pointer", INT2),  # 36-37
    ("vel_pointer", INT2),  # 38-39
    ("width_pointer", INT2),  # 40-41
    ("doppler_resolution", CODE2),  # 42-43
    ("vcp", INT2),  # 44-45
    ("spare_1", "8s"),  # 46-53
    ("spare_2", "2s"),  # 54-55
    ("spare_3", "2s"),  # 56-57
    ("spare_4", "2s"),  # 58-59
    ("nyquist_vel", SINT2),  # 60-61
    ("atmos_attenuation", SINT2),  # 62-63
    ("threshold", SINT2),  # 64-65
    ("spot_blank_status", INT2),  # 66-67
    ("spare_5", "32s"),  # 68-99
    # 100+  reflectivity, velocity and/or spectral width data, CODE1
)

# Table XI Volume Coverage Pattern Data (Message Types 5 & 7)
# pages 3-45 to 3-49 of ICD 2620002AA (halfword locations in comments).
# The former ("spare", "10s") is split into reserved halfwords and the VCP
# sequencing / supplemental halfwords. Total size is unchanged (22 bytes).
MSG_5 = (
    ("msg_size", INT2),  # HW 1
    ("pattern_type", CODE2),  # HW 2
    ("pattern_number", INT2),  # HW 3
    ("num_cuts", INT2),  # HW 4, 1 to 32
    # HW 5: VCP version in the upper byte, clutter map group in the lower byte
    ("clutter_map_group", INT2),
    ("doppler_vel_res", CODE1),  # HW 6 upper: 2 = 0.5 m/s, 4 = 1.0 m/s
    ("pulse_width", CODE1),  # HW 6 lower: 2 = short, 4 = long
    ("spare_7_8", "4s"),  # HW 7-8: reserved
    ("vcp_sequencing", CODE2),  # HW 9: RPG use, see _decode_vcp_sequencing
    ("vcp_supplemental", CODE2),  # HW 10: SAILS/MRLE/MPDA/base tilt flags
    ("spare_11", "2s"),  # HW 11: reserved
)

# Per-cut data, 23 halfwords (E1-E23) repeated for every elevation cut
MSG_5_ELEV = (
    ("elevation_angle", CODE2),  # E1, scaled by 360/65536 for value in degrees.
    ("channel_config", CODE1),  # E2 upper: 0 constant, 1 random, 2 SZ-2 phase
    ("waveform_type", CODE1),  # E2 lower: see _WAVEFORM_TYPES
    ("super_resolution", CODE1),  # E3 upper: bit field (bits 0-3)
    ("prf_number", INT1),  # E3 lower
    ("prf_pulse_count", INT2),  # E4
    ("azimuth_rate", CODE2),  # E5
    ("ref_thresh", SINT2),  # E6
    ("vel_thresh", SINT2),  # E7
    ("sw_thresh", SINT2),  # E8
    ("zdr_thres", SINT2),  # E9
    ("phi_thres", SINT2),  # E10
    ("rho_thres", SINT2),  # E11
    ("edge_angle_1", CODE2),  # E12
    ("dop_prf_num_1", INT2),  # E13
    ("dop_prf_pulse_count_1", INT2),  # E14
    # E15, formerly ("spare_1", "2s"). Same size, now per-cut SAILS/MRLE flags,
    # see _decode_elevation_supplemental.
    ("supplemental_data", CODE2),
    ("edge_angle_2", CODE2),  # E16
    ("dop_prf_num_2", INT2),  # E17
    ("dop_prf_pulse_count_2", INT2),  # E18
    ("spare_2", "2s"),  # E19 (EBC angle)
    ("edge_angle_3", CODE2),  # E20
    ("dop_prf_num_3", INT2),  # E21
    ("dop_prf_pulse_count_3", INT2),  # E22
    ("spare_3", "2s"),  # E23 (reserved)
)

# Table XVII-B Data Block (Descriptor of Generic Data Moment Type)
# pages 3-84 and 3-85
GENERIC_DATA_BLOCK = (
    ("block_type", "1s"),
    ("data_name", "3s"),  # VEL, REF, SW, RHO, PHI, ZDR
    ("reserved", INT4),
    ("ngates", INT2),
    ("first_gate", SINT2),
    ("gate_spacing", SINT2),
    ("thresh", SINT2),
    ("snr_thres", SINT2),
    ("flags", CODE1),
    ("word_size", INT1),
    ("scale", REAL4),
    ("offset", REAL4),
    # then data
)

# Table XVII-E Data Block (Volume Data Constant Type)
# page 3-86
VOLUME_DATA_BLOCK = (
    ("block_type", "1s"),
    ("data_name", "3s"),
    ("lrtup", INT2),
    ("version_major", INT1),
    ("version_minor", INT1),
    ("lat", REAL4),
    ("lon", REAL4),
    ("height", SINT2),
    ("feedhorn_height", INT2),
    ("refl_calib", REAL4),
    ("power_h", REAL4),
    ("power_v", REAL4),
    ("diff_refl_calib", REAL4),
    ("init_phase", REAL4),
    ("vcp", INT2),
    ("spare", "2s"),
)

# Table XVII-F Data Block (Elevation Data Constant Type)
# page 3-87
ELEVATION_DATA_BLOCK = (
    ("block_type", "1s"),
    ("data_name", "3s"),
    ("lrtup", INT2),
    ("atmos", SINT2),
    ("refl_calib", REAL4),
)

# Table XVII-H Data Block (Radial Data Constant Type)
# page 3-87
RADIAL_DATA_BLOCK = (
    ("block_type", "1s"),
    ("data_name", "3s"),
    ("lrtup", INT2),
    ("unambig_range", SINT2),
    ("noise_h", REAL4),
    ("noise_v", REAL4),
    ("nyquist_vel", SINT2),
    ("spare", "2s"),
)
