===========================================================================================
How to generate source code for the GOTM extension for Python using F2PY
===========================================================================================

-------------------------------------------------------------------------------------------

1. generate signature file gotm.pyf from GOTM's gotm.F90, time.F90 and our custom gui_util.f90

f2py.py -h gotm.pyf ../../src/gotm/gotm.F90 only: init_gotm time_loop clean_up : ../../src/util/time.F90 only: minn maxn : ../gui_util.f90 only: redirectoutput resetoutput getversion : -m gotm

(on Windows this is combined in f2py1.bat)

-------------------------------------------------------------------------------------------

2. the auto-generated gotm.pyf contains many unneeded module-level variables. We need only:

  time/minn
  time/maxn

Solve by removing all other variables from gotm.pyf. Also, use statements are obsolete and might be removed.

-------------------------------------------------------------------------------------------

3. in gotm.pyf, f2py currently chooses the wrong scope for time/minn, time/maxn; gotm.pyf lists 'private,public', but we need 'public'. Solve by removing 'private,' in gotm.pyf.

-------------------------------------------------------------------------------------------

4. generate C code (gotmmodule.c) and FORTRAN wrapper code (gotm-f2pywrappers2.f90) for Python extension module:

f2py.py gotm.pyf

Rename *.f90 to *.F90 (this is GOTM convention, which allows us to use GOTM's implicit makefile rules operating on *.F90). Rename gotmmodule.c to gotmmodule.cpp; it is now still C, but will be C++ when we have inserted our exception-handling code.

(on Windows this is combined in f2py4.bat)

-------------------------------------------------------------------------------------------

5. before we enter 'blocking' (i.e. time-consuming) FORTRAN code from Python, we must release the python Global Interpreter Lock (GIL). Releasing the GIL allows other Python threads to execute while the FORTRAN routine is running. To implement this:

In the auto-generated gotmmodule.cpp, look-up calls to FORTRAN code; these typically look like (*f2py_func)();

If this FORTRAN call is deemed expensive (in particular gotm/time_loop!), insert the Python macro Py_BEGIN_ALLOW_THREADS just above the call, and the macro Py_END_ALLOW_THREADS just below the call, each on a separate line on its own.

-------------------------------------------------------------------------------------------

6. *OBSOLETE* Now NetCDF is built so that the IFORT calling convention is supported by default; we do not need to change the default calling convention in the Visual Fortran project anymore.

In gotm-f2pywrappers2.f90, add at the start of every subroutine, i.e. below every 'subroutine somename(...)', add the lines

#ifdef WINDOWS
      !DEC$ ATTRIBUTES DEFAULT :: somename
#endif

where somename equals the name of the subroutine. This tells the Intel Fortran compiler to disable any non-default calling conventions (the current Windows project uses by default Compaq Visual Fortran compatibility, which wouold break our C-Fortran link without the above reset)

-------------------------------------------------------------------------------------------

7. *OBSOLETE* Now the Windows NetCDF lib is built so that the IFORT calling convention is supported by default; we do not need to change the default calling convention in the Windows Visual Fortran project anymore. The inerface block below still might be used, but the #ifdef block in any case must be removed.

In gotm-f2pywrappers2.f90, replace every line with 'external f2pysetupfunc' with an interface block of the like:

      interface
        subroutine f2pysetupfunc(...)
#ifdef WINDOWS
          !DEC$ ATTRIBUTES DEFAULT :: f2pysetupfunc
#endif
          ...
        end subroutine f2pysetupfunc
      end interface

For the first ..., you need to fill in the names of the arguments that will be given to the (C++) routine. For the second ..., the arguments must be declared (with type, often 'external') and the previously used name.

For an explanantion of the preprocessor (#...) lines, see item (6).

-------------------------------------------------------------------------------------------

8. *OBSOLETE* We now use the most recent NumPy on all platforms.

In gotmmodule.cpp, replace the import_array() line near the end of the file with

#ifdef OLD_NUMPY
  import_array();
#else
  import_array1(-1);
#endif

This is needed because the import_array syntax has changed with numpy versions. Currently the Windows build uses import_array(); (i.e. the old numpy), whereas the Linux build uses the new import_array1(-1);.

-------------------------------------------------------------------------------------------

9. In gotmmodule.cpp, add a function callable from FFortran that throws a C++ exception. E.g.:

void for_stop_core(void) {
	throw "GOTM library stopped\n";
}

All Fortran code should be changed to call this C++ function when a fatal error occurs, rather than 'stop'. The above name 'for_stop_core' is in fact the function called by IFORT when stop statements are encountered. Thus, changing the fortran code is not needed *as long as ifort is used as compiler*!

-------------------------------------------------------------------------------------------

10. Insert C++ exception handling code around calls to Fortran that could call the stop function of (9), in gotmmodule.cpp. Calls to Fortran from C++ look like

Py_BEGIN_ALLOW_THREADS
(*f2py_func)();
Py_END_ALLOW_THREADS

[Py_BEGIN_ALLOW_THREADS and Py_END_ALLOW_THREADS macros are optional, see step 5]

Change these to:

const char * strException = 0;
Py_BEGIN_ALLOW_THREADS
try {
        (*f2py_func)();
}
catch (const char * ex) {
	strException = ex;
}
Py_END_ALLOW_THREADS
if (strException!=0) PyErr_SetString(PyExc_Exception, strException);

11.

below lines

static FortranDataDef ...

replace all (void *) by (f2py_init_func)
