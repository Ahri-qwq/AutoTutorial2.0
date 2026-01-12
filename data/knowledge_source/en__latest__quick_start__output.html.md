# https://abacus.deepmodeling.com/en/latest/quick_start/output.html

> Source: https://abacus.deepmodeling.com/en/latest/quick_start/output.html

# Brief Introduction of the Output Files[#](#brief-introduction-of-the-output-files)

`OUT.suffix`

(Default one is `OUT.ABACUS`

). Here we give some simple descriptions.

*INPUT*[#](#input)

Different from `INPUT`

given by the users, `OUT.suffix/INPUT`

contains all parameters in ABACUS.


Note:`OUT.suffix/INPUT`

contains theactual parameters used in the calculation, including:


User-specified parameters(explicitly defined in your input file or command-line arguments, overriding default parameters).

System default parameters(automatically applied when not explicitly provided by the user).

For a complete list of input parameters, please consult this [instruction](../advanced/input_files/input-main.html).

*running_scf.log*[#](#running-scf-log)

`running_scf.log`

contains information on nearly all function calls made during the execution of ABACUS.

[KPT.info](http://KPT.info)[#](#kpt-info)

[KPT.info](http://KPT.info)

*eig.txt*[#](#eig-txt)

[istate.info](http://istate.info)’

*STRU.cif*[#](#stru-cif)

`.cif`

format structure file based on the input file `STRU`

, facilitating users to visualize with commonly used software.

*warning.log*[#](#warning-log)

The file contains all the warning messages generated during the ABACUS run.
