function findKthSmallestNumberRecursive(A,B,left,right,k)
     if left > right
           return -1
     end

     i = floor(Int32,(left+right)/2)
     j = k - i

     p = size(A,1)
     q = size(B,1)

     if j<0 || (j<=q && B[j+1] < A[i])
           return findKthSmallestNumberRecursive(A,B,left,i-1,k)
     elseif j>q || (j>0 && A[i] < B[j])
           return findKthSmallestNumberRecursive(A,B,i+1,right,k)
     else
           return i
     end
end

function findKthSmallestNumber(A,B,k)

      pos = findKthSmallestNumberRecursive(A,B,1,p,k)
      if pos>0
            return A[pos]
      end
      pos = findKthSmallestNumberRecursive(B,A,1,q,k)

      return B[pos]

end

arr1 = Array{Int32,1}(undef,10)
for i in 1:10
      arr1[i]=i
end

arr2 = Array{Int32,1}(undef,4)
for i in 1:4
      arr2[i]=i+3
end

print(findKthSmallestNumber(arr1,arr2,1))
