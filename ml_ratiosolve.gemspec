# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'ml_ratiosolve/version'

Gem::Specification.new do |spec|
  spec.name          = "ml_ratiosolve"
  spec.version       = MlRatiosolve::VERSION
  spec.authors       = ["Colin J. Fuller"]
  spec.email         = ["cjfuller@gmail.com"]
  spec.description   = %q{Methods for using maximum likelihood calculations to estimate parmeters of ratios of gaussian variates}
  spec.summary       = %q{Methods for using maximum likelihood calculations to estimate parmeters of ratios of gaussian variates}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files`.split($/)
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  #spec.add_dependency "nmatrix", ">= 0.0.9" #uncomment when required changes are in a gem release
  spec.add_dependency "trollop", "~> 2.0"
  spec.add_dependency "RubyInline"
  spec.add_development_dependency "bundler", "~> 1.3"
  spec.add_development_dependency "rake"
  spec.add_development_dependency "rspec"
end
