#! /usr/bin/env python
"""
autoEdit - A Python tool to automatically edit a set of files
           according to the specified user rules:
G. Kenway
"""

# Import modules
import os
import sys
import re

# Specify file extension
EXT = "_d.f90"

DIR_ORI = sys.argv[1]
DIR_MOD = sys.argv[2]

# Specify the list of LINE ID's to find, what to replace and with what
patt_modules = re.compile(r"(\s*use\s*\w*)(_d)\s*")
patt_module = re.compile(r"\s*module\s\w*")
patt_module_start = re.compile("(\s*module\s)(\w*)(_d)\s*")
patt_module_end = re.compile("(\s*end module\s)(\w*)(_d)\s*")
patt_subroutine = re.compile(r"\s*subroutine\s\w*")
patt_function = re.compile(r"\s*function\s\w*")

patt_subend = re.compile(r"\s*end\s*subroutine")
patt_funcend = re.compile(r"\s*end\s*function\n")

print("Directory of input source files  :", DIR_ORI)
print("Directory of output source files :", DIR_MOD)

useful_modules = [
    "bcroutines_d",
    "turbbcroutines_d",
    "utils_d",
    "flowutils_d",
    "walldistance_d",
    "bcpointers_d",
    "initializeflow_d",
    "turbutils_d",
    "sa_d",
    "fluxes_d",
    "solverutils_d",
    "residuals_d",
    "surfaceintegrations_d",
]
modSubToKeep = []

TDIA3X3_SOURCE = """

    subroutine tdia3x3(nb, ne, l, c, u, r)
        use constants
        implicit none
        integer(kind=inttype), intent(in) :: nb, ne
        real(kind=realtype), dimension(3, nb:ne), intent(inout) :: l, u, r
        real(kind=realtype), dimension(3, 3, nb:ne), intent(inout) :: c
        integer(kind=inttype) :: n
        real(kind=realtype) :: deti, r1, r2
        real(kind=realtype) :: ci11, ci12, ci13, ci21, ci22, ci23, ci31, ci32, ci33
        real(kind=realtype) :: f11, f12, f13, f21, f22, f23, f31, f32, f33

        do n = ne - 1, nb, -1
            ci11 = c(2, 2, n + 1) * c(3, 3, n + 1) - c(2, 3, n + 1) * c(3, 2, n + 1)
            ci12 = c(1, 3, n + 1) * c(3, 2, n + 1) - c(1, 2, n + 1) * c(3, 3, n + 1)
            ci13 = c(1, 2, n + 1) * c(2, 3, n + 1) - c(1, 3, n + 1) * c(2, 2, n + 1)
            ci21 = c(2, 3, n + 1) * c(3, 1, n + 1) - c(2, 1, n + 1) * c(3, 3, n + 1)
            ci22 = c(1, 1, n + 1) * c(3, 3, n + 1) - c(1, 3, n + 1) * c(3, 1, n + 1)
            ci23 = c(1, 3, n + 1) * c(2, 1, n + 1) - c(1, 1, n + 1) * c(2, 3, n + 1)
            ci31 = c(2, 1, n + 1) * c(3, 2, n + 1) - c(2, 2, n + 1) * c(3, 1, n + 1)
            ci32 = c(1, 2, n + 1) * c(3, 1, n + 1) - c(1, 1, n + 1) * c(3, 2, n + 1)
            ci33 = c(1, 1, n + 1) * c(2, 2, n + 1) - c(1, 2, n + 1) * c(2, 1, n + 1)
            deti = one / (c(1, 1, n + 1) * ci11 + c(1, 2, n + 1) * ci21 + c(1, 3, n + 1) * ci31)

            f11 = u(1, n) * ci11 * deti
            f12 = u(1, n) * ci12 * deti
            f13 = u(1, n) * ci13 * deti
            f21 = u(2, n) * ci21 * deti
            f22 = u(2, n) * ci22 * deti
            f23 = u(2, n) * ci23 * deti
            f31 = u(3, n) * ci31 * deti
            f32 = u(3, n) * ci32 * deti
            f33 = u(3, n) * ci33 * deti

            c(1, 1, n) = c(1, 1, n) - f11 * l(1, n + 1)
            c(1, 2, n) = c(1, 2, n) - f12 * l(2, n + 1)
            c(1, 3, n) = c(1, 3, n) - f13 * l(3, n + 1)
            c(2, 1, n) = c(2, 1, n) - f21 * l(1, n + 1)
            c(2, 2, n) = c(2, 2, n) - f22 * l(2, n + 1)
            c(2, 3, n) = c(2, 3, n) - f23 * l(3, n + 1)
            c(3, 1, n) = c(3, 1, n) - f31 * l(1, n + 1)
            c(3, 2, n) = c(3, 2, n) - f32 * l(2, n + 1)
            c(3, 3, n) = c(3, 3, n) - f33 * l(3, n + 1)

            r(1, n) = r(1, n) - f11 * r(1, n + 1) - f12 * r(2, n + 1) - f13 * r(3, n + 1)
            r(2, n) = r(2, n) - f21 * r(1, n + 1) - f22 * r(2, n + 1) - f23 * r(3, n + 1)
            r(3, n) = r(3, n) - f31 * r(1, n + 1) - f32 * r(2, n + 1) - f33 * r(3, n + 1)
        end do

        ci11 = c(2, 2, nb) * c(3, 3, nb) - c(2, 3, nb) * c(3, 2, nb)
        ci12 = c(1, 3, nb) * c(3, 2, nb) - c(1, 2, nb) * c(3, 3, nb)
        ci13 = c(1, 2, nb) * c(2, 3, nb) - c(1, 3, nb) * c(2, 2, nb)
        ci21 = c(2, 3, nb) * c(3, 1, nb) - c(2, 1, nb) * c(3, 3, nb)
        ci22 = c(1, 1, nb) * c(3, 3, nb) - c(1, 3, nb) * c(3, 1, nb)
        ci23 = c(1, 3, nb) * c(2, 1, nb) - c(1, 1, nb) * c(2, 3, nb)
        ci31 = c(2, 1, nb) * c(3, 2, nb) - c(2, 2, nb) * c(3, 1, nb)
        ci32 = c(1, 2, nb) * c(3, 1, nb) - c(1, 1, nb) * c(3, 2, nb)
        ci33 = c(1, 1, nb) * c(2, 2, nb) - c(1, 2, nb) * c(2, 1, nb)
        deti = one / (c(1, 1, nb) * ci11 + c(1, 2, nb) * ci21 + c(1, 3, nb) * ci31)
        r1 = r(1, nb)
        r2 = r(2, nb)
        r(1, nb) = deti * (ci11 * r1 + ci12 * r2 + ci13 * r(3, nb))
        r(2, nb) = deti * (ci21 * r1 + ci22 * r2 + ci23 * r(3, nb))
        r(3, nb) = deti * (ci31 * r1 + ci32 * r2 + ci33 * r(3, nb))

        do n = nb + 1, ne
            r(1, n) = r(1, n) - l(1, n) * r(1, n - 1)
            r(2, n) = r(2, n) - l(2, n) * r(2, n - 1)
            r(3, n) = r(3, n) - l(3, n) * r(3, n - 1)

            ci11 = c(2, 2, n) * c(3, 3, n) - c(2, 3, n) * c(3, 2, n)
            ci12 = c(1, 3, n) * c(3, 2, n) - c(1, 2, n) * c(3, 3, n)
            ci13 = c(1, 2, n) * c(2, 3, n) - c(1, 3, n) * c(2, 2, n)
            ci21 = c(2, 3, n) * c(3, 1, n) - c(2, 1, n) * c(3, 3, n)
            ci22 = c(1, 1, n) * c(3, 3, n) - c(1, 3, n) * c(3, 1, n)
            ci23 = c(1, 3, n) * c(2, 1, n) - c(1, 1, n) * c(2, 3, n)
            ci31 = c(2, 1, n) * c(3, 2, n) - c(2, 2, n) * c(3, 1, n)
            ci32 = c(1, 2, n) * c(3, 1, n) - c(1, 1, n) * c(3, 2, n)
            ci33 = c(1, 1, n) * c(2, 2, n) - c(1, 2, n) * c(2, 1, n)
            deti = one / (c(1, 1, n) * ci11 + c(1, 2, n) * ci21 + c(1, 3, n) * ci31)
            r1 = r(1, n)
            r2 = r(2, n)
            r(1, n) = deti * (ci11 * r1 + ci12 * r2 + ci13 * r(3, n))
            r(2, n) = deti * (ci21 * r1 + ci22 * r2 + ci23 * r(3, n))
            r(3, n) = deti * (ci31 * r1 + ci32 * r2 + ci33 * r(3, n))
        end do
    end subroutine tdia3x3
"""

for f in os.listdir(DIR_ORI):
    if f.endswith(EXT):
        # open original file in read mode
        file_object_ori = open(os.path.join(DIR_ORI, f), "r")
        print("\nParsing input file", file_object_ori.name)

        # read to whole file to string and reposition the pointer
        # at the first byte for future reading
        all_src = file_object_ori.read()
        file_object_ori.seek(0)

        # First we want to determine if it is a 'useful' module or a
        # 'useless' module. A useful module is one that has
        # subroutines in it.
        isModule = False
        hasSubroutine = False
        for line in file_object_ori:
            line = line.lower()
            if patt_module.match(line):
                isModule = True
            if patt_subroutine.match(line):
                hasSubroutine = True
            if patt_function.match(line):
                hasSubroutine = True

        # If we have a module, close the input and cycle to next file.
        if isModule and not hasSubroutine:
            file_object_ori.close()
            continue
        elif isModule and hasSubroutine:
            f = f.replace("_d", "_d")

        # open modified file in write mode
        file_object_mod = open(os.path.join(DIR_MOD, f), "w")

        # Go back to the beginning
        file_object_ori.seek(0)
        inSubroutine = False

        for line in file_object_ori:
            # Just deal with lower case string
            line = line.lower()

            # Replace _cb on calls
            if "_cd" in line:
                line = line.replace("_cd", "")

            if line.strip() == "external tdia3x3":
                line = ""

            # Replace _d modules with normal -- except for the useful
            # ones.
            m = patt_modules.match(line)
            if m:
                found = False
                for m in useful_modules:
                    if m in line:
                        found = True
                if found:
                    line = line.replace("_d", "_d", 1)
                else:
                    line = line.replace("_d", "")

            # # See if we need to modify the line with changing the
            # # module names
            # m = patt_module_start.match(line)
            # if m:
            #     line = 'module %s_d2\n'%m.group(2)

            # m = patt_module_end.match(line)
            # if m:
            #     line = 'end module %s_d2\n'%m.group(2)

            # Tapenade misses one function in inviscidupwindflux_d and we need to add it manually
            if patt_subroutine.match(line) and "inviscidupwindflux_d" in line:
                inSubroutine = True

            # If within the subroutine we just search for a very specific string append
            if inSubroutine and "use flowutils_d, only : etot" in line:
                line = line.strip("\n") + ", etot_d\n"

            if patt_subend.match(line):
                inSubroutine = False

            if f.lower() == "turbutils_d.f90" and line.strip() == "end module turbutils_d":
                file_object_mod.write(TDIA3X3_SOURCE)

            file_object_mod.write(line)

        # close the files
        file_object_ori.close()
        file_object_mod.close()

        # success message
        print(" Modified file saved", file_object_mod.name)
