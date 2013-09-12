#--
# error_bootstrapping_spec.rb
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

describe ErrorBootstrapping do 

  it "should generate random variables with the correct mean" do 
    ErrorBootstrapping.randnorm(10, 0).should eq 10.0
  end

  it "should generate random variables with reasonable statitstics" do
    vars = Array.new(100) { ErrorBootstrapping.randnorm(0, 1) }
    N[*vars].mean.to_f.abs.should < 0.4
    (N[*vars].std.to_f - 1).abs.should < 0.4
  end

  context "data generation" do
    it "should generate a data matrix of the correct size" do
      parameters = {mu: N[0,0], sig2: N[1,1], gamma: N[1]}
      x = N[[0.0], [1.0]]
      gd = ErrorBootstrapping.gen_data(parameters, x)
      gd.shape.should eq [2,1]
    end
  end
end

