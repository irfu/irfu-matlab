# Contributing to irfu-matlab

:+1::tada: Thank you kindly for taking the time to contribute, it is much appreciated! :tada::+1:


The following is a set of guidelines for contributing to [irfu-matlab](https://github.com/irfu/irfu-matlab/), which is hosted on GitHub.com

These are mostly guidelines, in other words they are not rules or laws or government regulations (*myndighetsföreskrifter* in Swedish). Use your best judgment, and feel free to propose changes to this document in a pull request.
If you are about to make any pull requests please request a merge to our [devel](https://github.com/irfu/irfu-matlab/tree/devel) branch, and in time these changes will propagate to our master branch when we are ready for a new release of irfu-matlab.


Coding style
------------
1. According to [Matlab Style Guidelines](http://www.datatool.com/prod02.htm "MATLAB Programming Style Guidelines")
2. Help follows matlab style, see [irf_template.m](https://github.com/irfu/irfu-matlab/blob/master/irf/irf_template.m)
3. Matlab code files checked for possible problems, using [checkcode](https://www.mathworks.com/help/matlab/ref/checkcode.html)

Git
---
1. Workflow according to https://nvie.com/posts/a-successful-git-branching-model/
2. Routines in [master](https://github.com/irfu/irfu-matlab/tree/master) branch should do not less and not more than written in their help!
3. Our [devel](https://github.com/irfu/irfu-matlab/tree/devel) branch includes latest delivered development changes for the next release.
4. Semantic release numbering https://semver.org/
5. Write clear and useful git commit messages https://chris.beams.io/posts/git-commit/

Common data types
-----------------
| Data type      | Comment  |
| -------------- | -------- |
| DataObject     | The same as dataset in CAA, equal to cdf file representation in Matlab. |
| TSeries        | TSeries class object, with useful methods for time series data. |
| VariableStruct | Structure including all variable data in original format and metadata, <br> usually extracted from DataObject |
| VariableIrf    | Variable as a simple structure in a default format for irfu-matlab, where: <br> Variable.t - TimeArray <br> Variable.[data] - data matrix of dimension [t x dim1 x dim2]... (typical example energy spectrograms) <br> Variable.[unit,label,dimunit,dimlabel,dimvec] <br> Variable.[vec] - data matrix of dimension [t x ndim], (for example ndim=3 for field vector) <br> Variable.[abs] - absolute value of vector in case exists Variable.vec |
| VariableMat    | Matrix where first column is time and the other columns are data [DEPRECATED] |


Time
----
The following time types are to be used:
* EpochTT (prefered). GenericTimeArray class with epoch internaly given in TT2000
* EpochUnix
* CDF epoch TT2000 in nanoseconds (int64)
* UNIX epoch in seconds
* CDF epoch in microseconds
* CDF epoch16
