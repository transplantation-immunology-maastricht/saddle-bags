# This file is part of EMBL-HLA-Submission.
#
# EMBL-HLA-Submission is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# EMBL-HLA-Submission is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with EMBL-HLA-Submission. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0 

# This file contains specifications for packaging of the MinION Extractor GUI
# As a standalone executable.  This file is meant to be used with pyinstaller
# http://www.pyinstaller.org/


# -*- mode: python -*-

block_cipher = None


a = Analysis(['AlleleSubmissionEMBL.py'],
             binaries=None,
             datas=None,
             hiddenimports=['six', 'packaging', 'packaging.requirements', 'packaging.version', 'packaging.specifiers', 'Tkinter', 'tkFileDialog', 'Tkconstants'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['tkinter'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='AlleleSubmissionEMBLWindows',
          debug=False,
          strip=False,
          upx=True,
          console=True )
