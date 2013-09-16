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

require 'csv'

require 'nmatrix'
require 'inline'

module InlineHelpers
    # 
    # The normal probability distribution.
    # @param  x [Numeric] The point at which to calculate the probability density
    # @param  m [Numeric] The mean of the distribution
    # @param  s2 [Numeric] The variance of the distribution
    # 
    # @return [Float] The probability density at the specified point
    # def normpdf(x, m, s2)
    #   1.0/Math.sqrt(2*Math::PI*s2)*Math.exp(-(x-m)**2/(2.0*s2))
    # end
    inline do |builder|
      builder.include '<math.h>'
      builder.c "
        double normpdf(const double x, const double m, const double s2) {
          const double pi = 3.1415926535897932384626;
          return 1.0/sqrt(2*pi*s2)*exp(-(x-m)*(x-m)/(2.0*s2));
        }
      "
    end

end

# 
# Functions for ML estimates of distribution parameters for ratios of Gaussian
# measurements.
# 
# This module aims to deal with the case where one needs to average 
# measurements that are made relative to an internal control, such that each 
# measurement to be averaged is the ratio of two Gaussian-distributed 
# variables, which is not itself Gaussian-distributed.
# 
# Each independent experiment is assumed to consist of I treatments, and there
# are N total experiments.  Each measurement x_{i,n} is assumed to be a Gaussian
# random variable with mean mu_i/gamma_n and variance s2_i/gamma_n. Each mu and
# s2 parameter depends only on the treatment and not the experiment number,
# whereas the scale factors gamma depend only on the experiment and not the
# treatment.
# 
# This is a more general case of the ratio-to-internal-control problem, and the
# latter is simply the case where the gamma_n are each set to the value of the
# internal control for experiment n.
# 
# @author Colin J. Fuller
# 
module MLRatioSolve
  self.extend InlineHelpers
  class << self

    # # 
    # # The normal probability distribution.
    # # @param  x [Numeric] The point at which to calculate the probability density
    # # @param  m [Numeric] The mean of the distribution
    # # @param  s2 [Numeric] The variance of the distribution
    # # 
    # # @return [Float] The probability density at the specified point
    # def normpdf(x, m, s2)
    #   1.0/Math.sqrt(2*Math::PI*s2)*Math.exp(-(x-m)**2/(2.0*s2))
    # end

    # 
    # Sets up the array of indices to skip from a string.
    # @param  skip_str [String] A string containing treatment, experiment index 
    #   pairs.  Pairs should have a comma internally and a colon between 
    #   successive pairs.  E.g. '1,2;0,1'
    # 
    # @return [void]
    def set_skip_indices(skip_str)
      if skip_str.empty?
        @skip_indices = []
      else
        skip_pairs = skip_str.split(":")
        skip_pairs.map! { |e| e.split(',').map { |ee| ee.to_i } }
        @skip_indices = skip_pairs
      end
    end

    # 
    # Gets and returns the list of (treatment, experiment) pairs to skip due to
    # missing data
    # 
    # @return [Array] An array of two element arrays containing treatment index,
    #     experiment index ordered pairs
    def skip_indices
      @skip_indices
    end

    # 
    # The log likelihood function for the datapoints given a set of distribution
    # parameters.
    # @param  gamma [NMatrix] The per-experiment scale factors (vector with 
    #   N elements)
    # @param  x [NMatrix] The data points, an IxN matrix.
    # @param  mu [NMatrix] The mean parameters for each treatment (vector with I 
    #   elements)
    # @param  sig2 [NMatrix] The variance parameters for each treatment (vector 
    #   with I elements)
    # 
    # @return [Float] The log likelihood of the data given the parameters
    def log_l_fct(gamma, x, mu, sig2)
      ltot = 0
      i = mu.size
      n = gamma.size
      i.times do |ii|
        n.times do |nn|
          next if skip_indices.include? [ii,nn]
          ltot += Math.log(normpdf(x[ii,nn], mu[ii]/gamma[nn], sig2[ii]/gamma[nn]**2))
        end
      end
      ltot
    end


    #
    # Given the per-experiment scale factors, calculate the ML estimate of the 
    # mean parameters
    # @param  gamma [NMatrix] The per-experiment scale factors (vector with N 
    #   elements)
    # @param  x [NMatrix] The data points, an IxN matrix.
    # 
    # @return [NMatrix] The ML estimates of the mean parameters (vector with I 
    #   elements)
    def calculate_mu_i(gamma, x)
      n = gamma.size
      i = x.shape[0]
      mu = NMatrix.zeros([i,1])
      if skip_indices.empty? then
        mu = x.dot(gamma) / n
      else
        i.times do |ii|
          count = 0
          n.times do |nn|
            next if skip_indices.include? [ii,nn]
            mu[ii] += gamma[nn]*x[ii,nn]
            count += 1
          end
          mu[ii] /= count
        end
      end

      mu
    end


    # 
    # Given the per-experiment scale factors, calculate the ML estimate of the 
    # variance parameters
    # @param  gamma [NMatrix] The per-experiment scale factors (vector with N 
    #   elements)
    # @param  x [NMatrix] The data points, an IxN matrix.
    # @param  mu [NMatrix] The precalculated ML estimates of the mean parameters 
    #   for the provided gamma (vector with I elements)
    # 
    # @return [NMatrix] The ML estimates of the variance parameters (vector with 
    #   I elements)
    def calculate_sig2_i(gamma, x, mu)
      n = gamma.size
      i = mu.size
      sig2 = NMatrix.zeros([i,1])
      
      i.times do |ii|
        count = 0
        n.times do |nn|
          next if skip_indices.include? [ii,nn]
          sig2[ii] += (gamma[nn]*x[ii,nn] - mu[ii])**2
          count += 1
        end
        sig2[ii] /= count
      end

      sig2
    end


    # 
    # Given the mean and variance parameters, calculate the ML estimate of the 
    # experimental scale factors
    # @param  x [NMatrix] The data points, an IxN matrix.
    # @param  mu [NMatrix] The mean parameters (vector with I elements)
    # @param  sig2 [NMatrix] The variance parameters (vector with I elements)
    # 
    # @return [NMatrix] The ML estimates of the per-expeirment scale factors (  
    #   vector with N elements)
    def calculate_gamma_n(x, mu, sig2)
      n = x.shape[1]
      i = mu.size
      gamma = NMatrix.zeros([n,1])
      gamma[0] = 1

      1.upto(n-1) do |nn|
        xm_over_s2 = 0
        x2_over_s2 = 0
        i.times do |ii|
          next if skip_indices.include? [ii,nn]
          xm_over_s2 += x[ii,nn]*mu[ii]/sig2[ii]
          x2_over_s2 += x[ii,nn]**2/sig2[ii]
        end
        gamma[nn] = (xm_over_s2 + Math.sqrt(xm_over_s2**2 + 4*i*x2_over_s2))/(2*x2_over_s2)
      end

      gamma
    end


    # 
    # Nicely prints the supplied distribution parameters
    # @param  mu [NMatrix] the mean parameters (vector with I elements)
    # @param  sig2 [NMatrix] the variance parameters (vector with I elements)
    # @param  gamma [NMatrix] the per-experiment scale factors (vector with N 
    #   elements)
    # @param  x [NMatrix] the experimental data (IxN matrix)
    # @param  iternum [NMatrix] the index of the current iteration
    # 
    # @return [void]
    def print_parameters(mu, sig2, gamma, x, iternum)
      puts "="*10
      puts "At iteration number #{iternum}:"
      puts "mu = #{mu.to_s}"
      puts "sig = #{sig2.map{ |e| Math.sqrt e }.to_s}"
      puts "gamma = #{gamma.to_s}"
      puts "Log likelihood = #{log_l_fct(gamma, x, mu, sig2)}"
      nil
    end


    # 
    # Does one iteration of estimation, taking a starting gamma, calculating mu 
    # and s2, and then recalculating gamma from these
    # @param  gamma_in [NMatrix] the per-experiment scale factors (vector with N 
    #   elements)
    # @param  x [NMatrix] the experimental data (IxN matrix)
    # @param  iternum [NMatrix] the index of the current iteration
    # 
    # @return [NMatrix] the new estimates of the per-experiment scale factors (
    #   vector with N elements)
    def do_single_iteration(gamma_in, x, iternum)
      mu = calculate_mu_i(gamma_in, x)
      sig2 = calculate_sig2_i(gamma_in,x,mu)
      gamma_out = calculate_gamma_n(x, mu, sig2)
      #print_parameters(mu, sig2, gamma_out, x, iternum)
      gamma_out
    end

    # 
    # Do a series of iterations of estimation from a supplied starting estimate.
    # @param  n_iter [Numeric] The maximum number of iterations to perform.
    # @param  x [NMatrix] the experimental data (IxN matrix)
    # @param  gamma_start [NMatrix] the initial guess for the gamma experimental 
    #   scale factors (vector with N elements)
    # @param  tol=nil [Numeric] if non-nil, the iterations will terminate early 
    #   if the absolute change in the likelihood between two successive 
    #   iterations is less than this
    # 
    # @return [Hash] the results of the iteration, containing keys for mu, sig2,
    #   gamma, and l, which are the mean, variance, scale factors, and
    #   likelihood, respectively, after iteration completed
    def do_iters_with_start(n_iter, x, gamma_start, tol=nil)
      gamma = gamma_start
      last_L = -1.0*Float::MAX
      n_iter.times do |it|
        gamma = do_single_iteration(gamma, x, it)
        m = calculate_mu_i(gamma,x)
        l = log_l_fct(gamma, x, m, calculate_sig2_i(gamma, x, m))
        if tol and (l-last_L).abs < tol then
          break
        end
        last_L = l
      end
      {mu: calculate_mu_i(gamma, x), sig2: calculate_sig2_i(gamma, x, calculate_mu_i(gamma, x)),
        gamma: gamma, l: log_l_fct(gamma, x, calculate_mu_i(gamma,x), calculate_sig2_i(gamma, x, calculate_mu_i(gamma,x)))}
    end

    # 
    # Successively yields points on an n-dimensional grid.
    # @param  max_depth [Integer] The number of dimensions.
    # @param  n_per [Integer] The number of points on the grid in each dimension.
    # @param  min [Float] The minimum value of gridded points.
    # @param  max [Float] The maximum value of gridded points.
    # @param  curr [Array] The current points in progress of being built.  
    #   This should be an array of max_depth elements.  Contents will be  
    #   overwritten.
    # @param  i [Integer] The current dimension.  Supply 0 initially.
    # 
    # @return [void]
    def grid_recursive(max_depth, n_per, min, max, curr, i, &b)
      if i == max_depth then
        yield curr
        return
      end
      step = (max-min)/n_per
      min.step(max, step) do |val|
        curr[i] = val
        grid_recursive(max_depth, n_per, min, max, curr, i+1, &b)
      end
      nil
    end

    # 
    # Runs iterations of estimation from starting guesses arrayed on a grid.
    # @param  n_starts [Numeric] Total number of guesses to generate
    # @param  n_gen [Integer] Number of dimensions
    # @param  n_iter [Integer] Maximum number of iterations to run per guess
    # @param  x [NMatrix] The datapoints (IxN matrix)
    # @param  tol=nil [Numeric] if non-nil, the iterations will terminate early 
    #   if the absolute change in the likelihood between two successive 
    #   iterations is less than this
    # 
    # @return [Hash] The result with the maximum likelihood among those found (
    #   see #do_iters_with_start for format)
    def grid_multiple_iters(n_starts, n_gen, n_iter, x, tol=nil)
      n_per = (n_starts**(1.0/n_gen)).ceil
      min_range = 0.01
      max_range = 5.0

      best = {l: -1.0*Float::MAX}
      # counter = 0

      grid_recursive(n_gen, n_per, min_range, max_range, Array.new(n_gen, 0.0), 0) do |gamma|
        result = do_iters_with_start(n_iter, x, N.new([gamma.size + 1,1],([1.0].concat(gamma))), tol)
        if result[:l] > best[:l] then
          best = result
        end
        # counter += 1
        # if best[:mu] then
        #   puts "Best solution found so far:"
        #   print_parameters(best[:mu], best[:sig2], best[:gamma], x, counter)
        # end
      end

      best
    end


    # 
    # Estimates the SEM for each treatment after normalization.
    # @param  result [Hash] The results hash as output by #do_iters_with_start
    # @param  norm_index [Integer] the index of the treatment to which to normalize
    # 
    # @return [NMatrix] The SEM estimates for each treatment (vector with I components)
    def ml_sem_estimate(result, norm_index)
      ni_skip_count = 0
      n = result[:gamma].size
      skip_indices.each do |si|
        if si[0] == norm_index then
          ni_skip_count += 1
        end
      end   
      prop_errors = result[:sig2].map.with_index do |e,i| 
        skip_count = 0
        skip_indices.each do |si|
          if si[0] == i then
            skip_count += 1
          end
        end
        sem2_i = e/(n-skip_count-1)
        sem2_norm = result[:sig2][norm_index]/(n-ni_skip_count-1)   
        sem2_i*1.0/result[:mu][norm_index]**2 + sem2_norm*(result[:mu][i]/result[:mu][norm_index]**2)**2
      end
      prop_errors.map { |e| Math.sqrt(e) }
    end


    # 
    # Reads the datapoints from a csv file.  Data should be arranged with 
    # columns as independent experiments and rows as treatments.
    # @param  fn [String] The filename from which to read the data.
    # 
    # @return [NMatrix] an IxN nmatrix containing the experimental data
    def read_data_from_file(fn)
      N[*(CSV.read(fn).map{ |e| e.map{ |i| i.to_f } })]
    end
  end
end




