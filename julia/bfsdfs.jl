#############################################################################
##
module bfsdfs

export gcd
export Node, BFS, DFs, pr, tree_print
export nodes, root  # test data

##  Euclid's algorithm in as a one-liner
gcd(a, b) = b == 0 ? a : gcd(b, a % b)

##  a tree type
struct Node
  id
  next::Array{Node}
end

##  a tree
nodes = [Node(i, []) for i in 1:7]
for (i,k) in pairs([3,4,4,5,6,6,6])
  i == k || push!(nodes[k].next, nodes[i])
end
root = nodes[6]

##  BFS
function BFS(x, visit)
  Q = [x]
  for y in Q
    visit(y)
    append!(Q, y.next)
  end
end

##  DFS
function DFS(x, visit)
  visit(x)
  for z in x.next
    DFS(z, visit)
  end
end

##  a visitor
pr(x) = print(x.id, ", ")

##  print as a tree
function tree_print(x, indent = "", first = true)
    first || print("\n", indent)
    print("-", x.id);
    first = true
    for c in x.next
        tree_print(c, indent * "  ", first)
        first = false
    end
end

end # module
