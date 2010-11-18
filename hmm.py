import sys, numpy as np
RENORM=1e-17

def forward_viterbi(obs, states, start_p, trans_p, emit_p):
   T = {}
   for state in states:
       ##          prob.           V. path  V. prob.
       T[state] = (start_p[state], [state], start_p[state])
   for output in obs:
       U = {}
       for next_state in states:
           total = 0
           argmax = None
           valmax = 0
           for source_state in states:
               (prob, v_path, v_prob) = T[source_state]
               p = emit_p[source_state][output] * trans_p[source_state][next_state]
               prob *= p
               v_prob *= p
               total += prob
               if v_prob > valmax:
                   argmax = v_path + [next_state]
                   valmax = v_prob
           U[next_state] = (total, argmax, valmax)
       T = U
   ## apply sum/max to the final states:
   total = 0
   argmax = None
   valmax = 0
   for state in states:
       (prob, v_path, v_prob) = T[state]
       total += prob
       if v_prob > valmax:
           argmax = v_path
           valmax = v_prob
   return (total, argmax, valmax)


class hmm(object):
   """Keeps track of states and solves optimal paths."""
   
   def _isTransitionMatrix(self, a):
      """Checks that a and b are transition matrices by making sure
      columns add rows are probabilites and that the number of hidden states
      match."""
      p=a.sum(1)
      if np.any(np.abs(p-1)>0.01):
         print 'WARNING transition matrix is not normalized'
         sys.exit(1)
      else:
         a=a/p

   def __init__(self, aa, bb):
      """Given transition probabilities and emission probabilities"""
      #Assert that data is correctly formated
      aa=np.asarray(aa); bb=np.asarray(bb)
      if aa.shape[0] != bb.shape[0]:
         print 'WARNING number of states in transition matrix and emission matrix disagree'
         sys.exit(2)
      if aa.ndim==2: 
         self._isTransitionMatrix(aa)
         self._isTransitionMatrix(bb)
         self.M = aa.shape[0]  #number of hidden states
         self.K = bb.shape[1]  #number of possible emitted states
      elif aa.ndim==3:
         for i in range(aa.shape[0]):
            self._isTransitionMatrix(aa[i])
            self._isTransitionMatrix(bb[i])
         self.M = aa.shape[1]  #number of hidden states
         self.K = bb.shape[2]  #number of possible emitted states
      self.a = aa
      self.b = bb

   def forward_backward(self, obs):
      """Solves for posterior probability of each."""
      M=self.M; K=self.K; N=len(obs); a=self.a; b=self.b; 
      alpha=np.zeros((N, M), np.float)
      beta =np.zeros((N, M), np.float)

      #Forward and backward step
      if a.ndim==3:
         alpha[0,:]=b[0,:,obs[0]]
      else:
         alpha[0,:]=b[:,obs[0]]
      beta[N-1,:]=1
      for tf in range(1,N):
         tb=N-1-tf
         if a.ndim==3:
            af=a[tf]; ab=a[tb]
            bf=b[tf]; bb=b[tb]
         else:
            af=a; ab=a
            bf=b; bb=b
         for j in range(M):
            for i in range(M): #sum step
               alpha[tf,j]+=alpha[tf-1,i]*af[i,j]*bf[j, obs[tf]]    #forward
               beta[tb,j] += beta[tb+1,i]*ab[j,i]*bb[i, obs[tb+1]]  #backward
         if alpha[tf,:].sum() < RENORM:
            alpha[tf,:]/=RENORM
         if beta[tb,:].sum() < RENORM:
            beta[tb,:]/=RENORM

      self.alpha=alpha
      self.beta=beta
      self.lhood=(alpha*beta).sum(1) #normalization factor
      self.pstate=((alpha*beta).T/self.lhood).T
      return self.pstate


if __name__ == '__main__':
    a=np.asarray([[0, .1, .2, .7, 0], 
                  [0,  0,  0,  1, 0],
                  [.1, 0, .8, .1, 0],
                  [.2, 0, .2, .4, .2],
                  [0,  0,  0, .3, .7]])
    
    b=np.asarray([[.2, 0, 0, .8, 0],
                  [0, 1., 0, 0, 0],
                  [.5, 0, 0, 0, .5],
                  [.2, 0, .6, .2, 0],
                  [.3, .2, .4, .1, 0]])
    y=np.asarray([3,3]*300)#,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3,1,3])
    probs=hmm(a,b)
    probs.forward_backward(y)
    print probs.pstate



