# Staden Package — Build Instructions (7 February 2011)

See [`NEWS`](NEWS) and [`ChangeLog`](ChangeLog) for version updates. Platform specific notes are preserved below.

## External Dependencies
To compile the Staden Package you will need the following third-party libraries. Some were historically maintained as separate projects (for example `io_lib`).

- `staden-io_lib` ≥ 1.13.8
- `tcl` ≥ 8.4
- `tk` ≥ 8.4
- `zlib` ≥ 1.2.x
- `liblzma` ≥ 4.999 *(optional – improves Gap5 compression)*
- `libpng` ≥ 1.2.x *(optional – required for Gap4 “Report Mutations” output)*
- `curl` ≥ 7.x *(required by io_lib)*
- `tklib` ≥ 0.3 *(tablelist widget)*
- `itcl` ≥ 3.2 *(optional – Gap4 auto-finishing GUI)*
- `itk` ≥ 3.2 *(optional – Gap4 auto-finishing GUI)*
- `iwidgets` ≥ 4.0 *(optional – Gap4 auto-finishing GUI)*

When your operating system separates runtime from development packages you must install both (for example `tcl8.5` **and** `tcl8.5-dev`).

### Other Toolchain Requirements
- C compiler (e.g. `gcc`, `icc`)
- C++ compiler (e.g. `g++`)
- FORTRAN compiler (e.g. `g77`). See `src/gap4/legacy_f2c.c` and `src/gap4/Makefile` if this is problematic.
- GNU `make`
- TexInfo *(documentation only)*
- Perl, Awk, Sed *(documentation only)*

> **Note**: When using `gfortran`, the GNU Autoconf system may miss the need for `-lgfortran` while linking. Add `LIBS=-lgfortran` to your `./configure` command if you observe link errors.

## Recommended Build Layout
The source tree expects out-of-tree builds.

```bash
cd src
mkdir build.myhost
cd build.myhost
../configure [options]
make
```

You may need to point `configure` at specific dependency install prefixes:

```
--with-io_lib=DIR
--with-zlib=DIR
--with-lzma=DIR
--with-libcurl=DIR
--with-tcl=DIR
--with-tk=DIR
--with-tclinclude=DIR
--with-tkinclude=DIR
--with-tklib=DIR
--prefix=DIR              # default: /usr/local/staden
```

Some platforms also require `LD_LIBRARY_PATH` (or the platform-specific equivalent) to include the runtime library directories during the build.

### Example Configure Summary
When configuration succeeds you should see a summary similar to:

```
External packages used:
curl:      via /usr/bin/curl-config
zlib:      DIR (system)
curses:    DIR (system)
liblzma:   DIR /software/badger/opt/xz_utils
io_lib:    via /software/badger/bin/io_lib-config
Tcl:       via /usr/lib/tcl8.4/tclConfig.sh
Tk:        via /usr/lib/tk8.4/tkConfig.sh
tklib:     /software/badger/opt/tcl_packages/tklib0.5
Iwidgets:  /usr/lib/iwidgets4.0.1
Itcl:      /usr/lib/itcl3.2
Itk:       /usr/lib/itk3.2
```

Do not be alarmed if optional packages (`iwidgets`, `itcl`, `itk`) are reported as `***NOT FOUND***`; they are only needed for the Gap4 auto-finishing GUI (`prefinish`).

## Platform Notes
### Linux / Debian
```bash
mkdir src/build.debian
cd src/build.debian
../configure \
    --with-tcl=/usr/lib/tcl8.4 --with-tk=/usr/lib/tk8.4 \
    --with-tklib=/software/badger/opt/tcl_packages/tklib0.5 \
    --with-lzma=/software/badger/opt/xz_utils \
    --with-io_lib=/software/badger \
    --prefix=$HOME/staden.install
```

### Linux / Ubuntu 10
Packages used (list may be incomplete):

```bash
apt-get install g++
apt-get install zlib1g-dev
apt-get install tk-dev
apt-get install liblzma-dev
apt-get install tklib
apt-get install libcurl3-gnutls-dev
apt-get install libncurses5-dev
```

### Linux / Fedora
```bash
yum install gcc
yum install gcc-c++
yum install zlib-devel
yum install tk-devel
yum install tklib
yum install curl-devel
yum install ncurses-devel
yum install libpng-devel
yum install libXt-devel  # for X11/Intrinsic.h
```

If SELinux prevents loading shared libraries (`cannot restore segment prot after reloc: Permission denied`), run:

```bash
chcon -t texrel_shlib_t /PREFIX/lib/staden/lib*
```

Replace `/PREFIX` with the value passed to `--prefix`. Alternatively, relax SELinux temporarily with `setenforce 0` (not generally recommended).

On 64-bit systems with 32-bit compatibility libraries, Fedora and similar distributions may prefer 32-bit builds. Override the search paths explicitly, e.g.

```bash
../configure --with-tcl=/usr/lib64 --with-tk=/usr/lib64
```

### Linux / CentOS
```bash
yum install xz-devel
yum install tk-devel
yum install zlib-devel
```

Other packages mirror those listed for Fedora. Refer to the SELinux notes above when mixing 32-bit and 64-bit libraries.

### Microsoft Windows
1. **Main source tree** – The same `configure` script can be used under MinGW/MSYS. Expect to install the GNU toolchain manually. You may also need to run `mingw-get install libz` to obtain headers (behaviour has varied over time).
2. **`src/windows/run`** – Build the small launcher that substitutes the Unix shell scripts (e.g. `gap4`, `trev`). It sets environment variables before invoking `tclsh` on the relevant GUI program.
3. **WiX installer** – Install WiX to produce `.msi` installers:
   ```bash
   cd src/windows/wix
   perl generate_wxs.pl c:/jkb/staden/install_root > staden_files.wxi
   DISTROOT=/c/jkb/staden/install_root candle staden.wxs
   light -ext WixUIExtension staden.wixobj
   ```
   Strawberry Perl is recommended because the MSYS Perl sometimes rewrites Windows paths incorrectly.

### macOS
#### Modern guidance (31 Oct 2013)
- Build Tcl/Tk (e.g. 8.6.1) from the `unix` directory, set `MACOSX_DEPLOYMENT_TARGET=10.5`, and configure with `CFLAGS=-arch i386 -arch x86_64`.
- After `make` and `make install`, copy `license.terms` into `$HOME/sys/lib/{tcl,tk}*` to satisfy redistribution terms.
- Build `tcllib` and `tklib` similarly, copying their licence files alongside the binaries.
- For `io_lib` add `--disable-dependency-tracking` and `--with-libcurl=$HOME/sys`.

#### Earlier guidance
- Tk comes in two forms: native Aqua frameworks and an X11 build. The Staden Package requires the X11 version.
- If `configure` finds Aqua first, point to the X11 installs manually, e.g.
  ```bash
  ./configure --with-tcl=/opt/local/lib --with-tk=/opt/local/lib
  ```
  (MacPorts layout shown.)
- Historical MacPorts package versions: `tcllib 1.11`, `tklib 0.5`, `tcl/tk 8.5.9`.

#### Building universal binaries
```bash
cd src/tcl8.5.9/unix
CFLAGS="-O -arch i386 -arch x86_64 -arch ppc -arch ppc64" ./configure --prefix=/nfs/sam_scratch/jkb/src/tcltk_inst
make && make install

cd ../../src/tk8.5.9/unix
CFLAGS="-O -arch i386 -arch x86_64 -arch ppc -arch ppc64" ./configure --prefix=/nfs/sam_scratch/jkb/src/tcltk_inst --with-tcl=/nfs/sam_scratch/jkb/src/tcltk_inst8.5/lib
make && make install
```

When compiling Staden remember to mirror `CFLAGS` into `CXXFLAGS` because some components are in C++.

#### Debugging universal binaries
- Inspect architectures with `otool` and `lipo`:
  ```bash
  otool -arch all -h libstaden-read.dylib
  otool -L libstaden-read.dylib
  lipo -info libstaden-read.dylib
  ```
- Print runtime linking details while launching Gap5:
  ```bash
  DYLD_PRINT_LIBRARIES=1 DYLD_PRINT_RPATHS=1 DYLD_PRINT_OPTS=1 ./bin/gap5
  ```

## Documentation and Course Notes
Historical course material moved into a separate archive `staden_doc-2.0.0b8-src.tar.gz` (and later releases). Download and consult that package for tutorials.

