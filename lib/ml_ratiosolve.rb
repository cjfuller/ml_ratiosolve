#--
# ml_ratiosolve.rb
# Copyright (c) 2013 Colin J. Fuller
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the Software), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED AS IS, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#++

require "ml_ratiosolve/version"
require 'ml_ratiosolve/ml_ratiosolve'
require 'ml_ratiosolve/error_bootstrapping'

module MLRatioSolveBin

  def self.go(opts)
    n_starts = opts[:n_starts]
    n_iter = opts[:n_iter]
    n_for_bootstrap = opts[:n_bootstrap]
    tol = opts[:tol]
    ci_level = opts[:ci_level]
    norm_index = opts[:norm_index]

    MLRatioSolve.set_skip_indices(opts[:skip])

    x = MLRatioSolve.read_data_from_file(opts[:file])

    n_gammas_to_fit = x.shape[1] - 1

    best = MLRatioSolve.grid_multiple_iters(n_starts, n_gammas_to_fit, n_iter, x, tol)

    puts "Best solution found: "
    puts "mu: #{best[:mu]/best[:mu][norm_index]}"
    puts "sig2: #{best[:sig2].map{ |e| Math.sqrt(e)/best[:mu][norm_index] }.to_s}"
    puts "gamma: #{best[:gamma]}"
    puts "log l: #{best[:l]}"

    puts "Error estimate:"
    puts MLRatioSolve.ml_sem_estimate(best, norm_index)

    #sim_results = ErrorBootstrapping.estimate_with_gen_data(n_for_bootstrap, best, x, n_iter, tol)
    #ci_lower, ci_upper = ErrorBootstrapping.bootstrap_ci(sim_results, ci_level)
    # puts "boostrapped #{ci_level*100}% confidence interval: "
    # puts ci_lower.to_a.join(", ")
    # puts ci_upper.to_a.join(", ")
  end

end
