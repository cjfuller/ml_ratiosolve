#--
# ml_ratiosolve_spec.rb
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

require 'ml_ratiosolve'

describe MLRatioSolve do

  before :each do
    MLRatioSolve.set_skip_indices ""
  end

  it "should calculate a normal PDF" do
    MLRatioSolve.normpdf(0, 0, 1).should be_within(1.0e-6).of 0.3989423
  end

  it "should set and retrieve skip indices" do
    MLRatioSolve.set_skip_indices("0,1:0,2")
    MLRatioSolve.skip_indices.should eq [[0,1],[0,2]]
  end

  it "should calculate the log likelihood" do
    gamma = N.new([2,1],[1.0, 0.5])
    mu = N.new([2,1], [0.0, 2.0])
    sig2 = N.new([2,1], [1.0, 16.0])
    x = N.new([2,2], [0.0, 0.0, 2.0, 2.0])
    expected = N[Math.log(0.199471140200716),
                Math.log(0.398942280401433),
                Math.log(0.0997355701003582),
                Math.log(0.0483335146003562)].sum.to_f
    MLRatioSolve.log_l_fct(gamma, x, mu, sig2).should be_within(1.0e-6).of expected
  
  end

  it "should calculate the ML mean estimate with uniform gamma" do
    x = N.new([1,3], [0,3,6], dtype: :float64)
    g = N.new([3,1], [1,1,1], dtype: :float64)
    result = MLRatioSolve.calculate_mu_i(g,x)
    result.shape.should eq [1,1]
    result[0].should eq 3.0
  end

  it "should calculate the ML mean estimate with nonuniform gamma" do
    x = N.new([1,2], [1,2], dtype: :float64)
    g = N.new([2,1], [4,1], dtype: :float64)
    result = MLRatioSolve.calculate_mu_i(g,x)
    result.shape.should eq [1,1]
    result[0].should eq 3.0
  end

  it "should calculate the ML variance estimate" do
    x = N.new([1,3], [0,3,6], dtype: :float64)
    g = N.new([3,1], [1,2,1], dtype: :float64)
    result = MLRatioSolve.calculate_sig2_i(g,x, MLRatioSolve.calculate_mu_i(g,x))
    result.shape.should eq [1,1]
    result[0].should eq N[0,6,6].variance.to_f*(2.0/3.0) #final fraction is to correct sample vs. population
  end

  it "should calculate a single gamma for fixed mean and standard deviation parameters" do
    mu = N.new([2,1], [0.0, 2.0])
    sig2 = N.new([2,1], [1.0, 16.0])
    x = N.new([2,2], [0.0, 0.0, 2.0, 2.0])
    MLRatioSolve.calculate_single_gamma(0, x, mu, sig2).should be_within(1.0e-6).of 3.3722813
  end

  it "should correctly normalize the set of gammas" do
    mu = N.new([2,1], [0.0, 2.0])
    sig2 = N.new([2,1], [1.0, 16.0])
    x = N.new([2,2], [0.0, 0.0, 2.0, 2.0])
    result = MLRatioSolve.calculate_gamma_n(x, mu, sig2)
    g0 = MLRatioSolve.calculate_single_gamma(0,x,mu,sig2)
    g1 = MLRatioSolve.calculate_single_gamma(1,x,mu,sig2)

    result.should eq N[[g0],[g1]]/g0
  end

  it "should perform a single iteration" do
    pending
  end

  it "should not regress on results of a set of iterations" do
    pending
  end

  context "grid search" do
    it "should grid the correct number of iterations" do
      pending    
    end

    it "should not regress on multiple gridded iterations" do
      pending
    end
  end

  context "low variance solution search" do
    it "should test the correct low variance solutions" do
      pending
    end

    it "should calculate the starting means for the low variance case" do
      pending
    end

    it "should calculate the starting variances for the low variance case" do
      pending
    end

    it "should find a permutation such that the 0th entry of a given low variance treatment is not skipped" do
      pending
    end

    it "should be able to invert a permutation matrix" do
      pending
    end

    it "should be able to permute rows of a matrix" do
      pending
    end

  end

  context "error assessment" do
    it "should correctly calculate the SEM estimate" do
      pending
    end
  end

end

