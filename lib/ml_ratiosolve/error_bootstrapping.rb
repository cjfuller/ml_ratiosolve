#--
# error_bootstrapping.rb
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

# 
# Methods for using parametric bootstrapping to estimate confidence intervals
# for the ML ratio estimation.
# 
# @author Colin J. Fuller
# 
module ErrorBootstrapping
  class << self

    # 
    # Generate a random normal variate using the Box-Muller transform
    # @param  mu [Numeric] The desired mean
    # @param  s2 [Numeric] The desired variance
    # 
    # @return [Float] A random variate from the normal distribution with 
    #   supplied parameters
    def randnorm(mu, s2)
      ru0 = rand
      ru1 = rand
      ([Math.sqrt(-2.0 * Math.log(ru0)) * Math.cos(2*Math::PI*ru1), Math.sqrt(-2.0 * Math.log(ru0)) * Math.sin(2*Math::PI*ru1)].sample)*Math.sqrt(s2) + mu
    end


    # 
    # Generate a set of simulated data consisting of random numbers drawn from 
    # the distributions with the supplied parameters
    # @param  parameters [Hash] A hash containing the mean, variance, and scale
    #   parameters formatted like the output from 
    #   MLRatioSolve::do_iters_with_start
    # @param  x [NMatrix] The original experimental data
    # 
    # @return [NMatrix] A matrix of simulated data with the same dimensions as 
    # x.  Any skipped data points (as returned by MLRatioSolve::skip_indices) 
    # will be set to 0 here.
    def gen_data(parameters, x)
      mu = parameters[:mu]
      sig2 = parameters[:sig2]
      gamma = parameters[:gamma]
      size_i = x.shape[0]
      size_n = x.shape[1]

      sim_data_out = N.zeros_like(x)

      size_i.times do |i|
        size_n.times do |n|
          next if MLRatioSolve.skip_indices.include?([i,n])
          sim_data_out[i,n] = randnorm(mu[i]/gamma[n], sig2[i]/gamma[n]**2)
        end
      end

      sim_data_out
    end


    # 
    # Re-estimate distribution parameters by generating simulated data a number 
    # of times and performing the iterative estimation in MLRatioSolve
    # @param  n_gen [Numeric] number of datasets to simulate
    # @param  parameters [Hash] A hash containing the parameters from the estimation run on the simulated data.
    # @param  x [NMatrix] The original experimental data
    # @param  n_iter [Numeric] max number of iterations per estimate
    # @param  tol=nil [Numeric] if non-nil, the iterations will terminate early 
    #   if the absolute change in the likelihood between two successive 
    #   iterations is less than this
    # 
    # @return [Array] An array containing n_gen hashes, each of which is the 
    #   result of the estimation run on one simulated dataset.
    def estimate_with_gen_data(n_gen, parameters, x, n_iter, tol=nil)
      result = Array.new(n_gen) { nil }
      n_gen.times do |i|
        xi = gen_data(parameters, x)
        result[i] = MLRatioSolve.do_iters_with_start(n_iter, xi, parameters[:gamma], tol)
      end
      result
    end


    # 
    # Calculate a bootstrapped confidence interval from the output of 
    # #estimate_with_gen_data.
    # @param  results [Array] An array of hashes as output from 
    #   #estimate_with_gen_data
    # @param  level [Float] A number between 0 and 1 that is the level of the 
    #   confidence interval.  For instance, a value of 0.95 will lead to 
    #   calculation of the bounds on the central 95% of values.
    # 
    # @return [Array] an array of two NMatrix objects, the lower bound on the   
    #   CI and the upper bound on the CI
    def bootstrap_ci(results, level)
      means = results.map { |e| e[:mu].to_a }
      means = N[*means]
      size_i = means.shape[1]
      size_sim = means.shape[0]
      ci_lower = N.new([size_i], 0.0)
      ci_upper = N.zeros_like(ci_lower)
      size_i.times do |i|
        means_i = means[0...size_sim, i].to_a
        means_i.flatten!
        means_i = means_i.select { |e| e.finite? }
        means_i.sort!
        lower_ci_index = ((1.0-level)/2.0 * means_i.length).ceil
        upper_ci_index = ((1.0 - (1.0-level)/2.0) * means_i.length).floor
        ci_lower[i] = means_i[lower_ci_index]
        ci_upper[i] = means_i[upper_ci_index]
      end
      [ci_lower, ci_upper]
    end

  end
end
