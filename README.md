# MlRatiosolve

A tool for parameter estimations for (not-quite-)ratio problems.  

In scientific data, frequently there is uncontrolled experiment-to-experiment variability such that quantitiative comparison between experiments is only meaningful when data are expressed as a ratio to an internal control that is part of each experiment.  However, even if the values of the control and other treatments would be Gaussian distributed, this ratio is not, making accurate estimates of means and errors difficult.

MLRatiosolve uses the fact that this problem contains more information than just the ratios to obtain better estimates of distribution parameters than is possible for a ratio problem.


## Installation

  Dependencies:

  MLRatioSolve requires the [NMatrix linear algebra library](https://github.com/SciRuby/nmatrix), which will be installed automatically.  NMatrix, however, requires that you first install ATLAS.  If you don't have this installed already, see the instructions at https://github.com/SciRuby/nmatrix/wiki/Installation

  `$ gem install ml_ratiosolve`

  Or, if you have ruby installed globally on your computer:

  `$ sudo gem install ml_ratiosolve`

## Usage

MLRatiosolve takes a CSV file of data as input.  Run `ml_ratiosolve --help` from the command line for further information and additional options.

## Contributing

1. Fork it
2. Create your feature branch (`git checkout -b my-new-feature`)
3. Commit your changes (`git commit -am 'Add some feature'`)
4. Push to the branch (`git push origin my-new-feature`)
5. Create new Pull Request
