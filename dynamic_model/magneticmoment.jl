using LinearAlgebra  

function f(i_1,i_2,i_3,B)
 #parameters   i_k:current , B:geomagnetism
 
 #constants
  μ=1#magnetic permeability
  n=1#the number of turns
  S=1#section area
  m_r=[1,1,1]#regidual magnetic moment

 #Calculate magnetic moments 
  m_1=μ*n*S*i_1*[1,0,0]
  m_2=μ*n*S*i_2*[0,1,0]
  m_3=μ*n*S*i_3*[0,0,1]

 #Calculate magnetic torques
  T_1=cross(m_1,B)
  T_2=cross(m_2,B)
  T_3=cross(m_3,B) 
  T_l=cross(m_r,B)

 #Find out all torgue
 T_1 + T_2 + T_3 + T_l

end