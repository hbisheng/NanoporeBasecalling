import numpy as np
import pickle

def sigmoid(x):
  return 1 / (1 + np.exp(-x))

class OutLayer:
  def __init__(self):
    self.n_params = 2
    self.params = [None, None]

  def calc(self, input):
    otmp = np.dot(input, self.params[0]) + self.params[1]
    e_x = np.exp(otmp - otmp.max(axis=1, keepdims=True))
    return e_x / e_x.sum(axis=1, keepdims=True)

class SimpleLayer:
  def __init__(self):
    self.n_params = 5
    self.params = [None for i in range(5)]
 
  def conv(self, input, w, b, width=7):
    output = np.zeros((input.shape[0], w.shape[0]))
    mid = width/2
    for i in range(width):
      wx = w[:,:,i,0].T
      oo = np.dot(input, wx)
      if i < mid:
        output[:-(mid-i)] += oo[mid-i:]
      elif i == mid:
        output += oo
      else:
        output[i-mid:] += oo[:-(i-mid)]

    output += b
    return output

  def calc(self, input):
    state = self.params[4]

    c1 = np.tanh(self.conv(input, self.params[0], self.params[1]))
    f1 = sigmoid(self.conv(input, self.params[2], self.params[3]))

    output = np.zeros((len(input), self.params[4].shape[0]), dtype=np.float32)
    for i in range(len(input)):
      state = f1[i] * state + (1 - f1[i]) * c1[i]
      output[i] = state
    return np.array(output)

class BiSimpleLayer:
  def __init__(self):
    self.fwd = SimpleLayer()
    self.bwd = SimpleLayer()

  def calc(self, input):
    return np.concatenate([self.fwd.calc(input), self.bwd.calc(input[::-1])[::-1]],
                          axis=1)

class Rnn:
  def __init__(self):
    pass

  def tester(self, input):
    input = input[0]
    l1 = self.layer1.calc(input)
    l2 = self.layer2.calc(l1)
    l3 = self.layer3.calc(l2)
    l4 = self.layer4.calc(l3)
    return [self.output1.calc(l4)], [self.output2.calc(l4)]

  def predict(self, input):
    l1 = self.layer1.calc(input)
    l2 = self.layer2.calc(l1)
    l3 = self.layer3.calc(l2)
    l4 = self.layer4.calc(l3)
    return self.output1.calc(l4), self.output2.calc(l4)

  def debug(self, input):
    l1 = self.layer1.calc(input)
    l2 = self.layer2.calc(l1)
    l3 = self.layer3.calc(l2)
    return l1, l2, l3

  def load(self, fn):
    with open(fn, "rb") as f:
      self.layer1 = BiSimpleLayer()
      for i in range(5):
        self.layer1.fwd.params[i] = pickle.load(f)
      for i in range(5):
        self.layer1.bwd.params[i] = pickle.load(f)
      self.layer2 = BiSimpleLayer()
      for i in range(5):
        self.layer2.fwd.params[i] = pickle.load(f)
      for i in range(5):
        self.layer2.bwd.params[i] = pickle.load(f)
      self.layer3 = BiSimpleLayer()
      for i in range(5):
        self.layer3.fwd.params[i] = pickle.load(f)
      for i in range(5):
        self.layer3.bwd.params[i] = pickle.load(f)
      self.layer4 = BiSimpleLayer()
      for i in range(5):
        self.layer4.fwd.params[i] = pickle.load(f)
      for i in range(5):
        self.layer4.bwd.params[i] = pickle.load(f)
      self.output1 = OutLayer()
      self.output2 = OutLayer()
      for i in range(2):
        self.output1.params[i] = pickle.load(f)
      for i in range(2):
        self.output2.params[i] = pickle.load(f)
