# Virtual Soil Simulator (VSS)

Richards' equation solver for OpenFOAM, incorporating the film-flow phenomenon in unsaturated hydraulic conductivity and soil water retention modelling.

---

## Requirements

- **OpenFOAM** (ESI/OpenCFD edition from [openfoam.com](https://www.openfoam.com/)).  Software tested under: OF2312, OF2406, OF2412 and OF2512.
- **Boost** development headers
- **Eigen3** development headers
- **Git**

---

## Installation
First install the OpenFOAM. Below procedure for two example distributions is provided.

### CentOS 9 Stream

```bash
# OpenFOAM
sudo dnf install -y dnf-plugins-core epel-release
sudo dnf copr enable -y openfoam/openfoam
sudo dnf install -y openfoam2512-default

# VSS dependencies
sudo dnf install -y git boost-devel eigen3-devel
```

### Ubuntu 26.04

```bash
# OpenFOAM
curl -s https://dl.openfoam.com/add-debian-repo.sh | sudo bash
sudo apt-get update
sudo apt-get install -y openfoam2512-default

# VSS dependencies
sudo apt-get install -y git libboost-all-dev libeigen3-dev
```

### Other Linux distributions

Follow the OpenFOAM installation guide for your platform:  
<https://gitlab.com/openfoam/core/openfoam/-/wikis/precompiled>

Then install Boost and Eigen3 development headers using your distribution's package manager.

---

## VSS installation
### Download the sources

```bash
git clone https://github.com/klamorski/vss.git
cd vss
```

---

### Compilation


#### 1. Configure the build

Default build options shall be appriopriate for most cases. If not, open `src/build.cfg` in a text editor. The key settings are:

| Variable | Values | Description |
|---|---|---|
| `installGlobal` | `true` / `false` | Where compiled binaries are installed (see below) |
| `buildMode` | `opt` / `debug` | Optimised (production) or debug build |
| `noCompileCPU` | integer | Number of CPU cores used for compilation |
| `makeType` | `libso` / `lib` | Shared (`.so`) or static (`.a`) library |

#### Local installation (default, no root required)

Set `installGlobal=false` in `build.cfg`. Binaries are placed in your personal OpenFOAM user directory (`$FOAM_USER_APPBIN` / `$FOAM_USER_LIBBIN`), typically `~/OpenFOAM/<user>-<version>/platforms/.../`.  
Only your user account will have access to the solver.

#### System-wide installation (all users, requires root)

Set `installGlobal=true` in `build.cfg`. Binaries are placed in the shared OpenFOAM directories (`$FOAM_APPBIN` / `$FOAM_LIBBIN`).  Write permissions to system folders where OpenFOAM is installed are needed in that case.

#### 2. Activate the OpenFOAM environment

Before building or running anything, the OpenFOAM environment must be initialised in your shell session. This makes all OpenFOAM tools (`wmake`, `wmakeLnInclude`, solvers, etc.) available on the `PATH`. If global istallation is used, issue `su` command before activating OpenFOAM environment.

```bash
openfoam2512        # This is for OpenFOAM v. 2512 call similar command loading 
                    # other OpenFOAM version if needed
```

> You need to repeat this every time you open a new terminal.

#### 3. Compile the software

```bash
cd vss/src
./Allwmake
```

#### 4. Cleaning a previous build

To remove compiled objects from source tree and restart from scratch:

```bash
cd vss/src
./Allwclean
```

---

## Running the solver

After a successful build, activate the OpenFOAM environment and call the solver from within a case directory:

```bash
openfoam2512
cd path/to/your/case
vssFoam
```

Example tutorial cases are provided in the `tutorial/` directory. Each case contains a `run.sh` script and a `clean.sh` script.
---

## License

VSS is distributed under the **GNU General Public License v3.0 or later** (GPL-3.0-or-later), in accordance with the license of OpenFOAM on which it depends. See the [LICENSE](LICENSE) file for the full text.

The following third-party libraries vendored in `src/library/external/` are distributed under their own permissive licenses, which are compatible with GPL v3:

| Library | Location | License |
|---|---|---|
| [fast-cpp-csv-parser](https://github.com/ben-strasser/fast-cpp-csv-parser) | `src/library/external/csv/` | BSD 3-Clause |
| [libInterpolate](https://github.com/CD3/libInterpolate) | `src/library/external/interpolate/` | MIT |
