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

  it "should calculat a normal PDF" do
    MLRatioSolve.normpdf(0, 0, 1).should be_within(0.01).of 0.4
  end

  it "should calculate the ML mean estimate with uniform gamma" do
    x = N.new([1,3], [0,3,6], :float64)
    g = N.new([3,1], [1,1,1], :float64)
    result = MLRatioSolve.calculate_mu_i(g,x)
    result.shape.should eq [1,1]
    result[0].should eq 3.0
  end

  it "should calculate the ML mean estimate with nonuniform gamma" do
    x = N.new([1,2], [1,2], :float64)
    g = N.new([2,1], [4,1], :float64)
    result = MLRatioSolve.calculate_mu_i(g,x)
    result.shape.should eq [1,1]
    result[0].should eq 3.0
  end

  it "should calculate the ML variance estimate" do
    x = N.new([1,3], [0,3,6], :float64)
    g = N.new([3,1], [1,2,1], :float64)
    result = MLRatioSolve.calculate_sig2_i(g,x, MLRatioSolve.calculate_mu_i(g,x))
    result.shape.should eq [1,1]
    result[0].should eq N[0,6,6].variance.to_f
  end

  it "should grid the correct number of iterations" do
    pending    
  end

end

