a = [1:5]
function kill_a()
  println(a)
  a[5] = 70000000
end
println(a)
kill_a()
println(a)

