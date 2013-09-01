--[[
plotting script

plot.lua <filename.txt> <dim> <field>
--]]

local filename, dim, field = ...
dim = assert(tonumber(dim), 'expected a valid dimension')

local useHistory = false
local useColor = true 

local f = assert(io.open(filename, 'rb'), "failed to open file")
local l = f:read('*l')
f:close()

assert(l:sub(1,1) == '#', "expected a comment on the first line")
l = l:sub(2)

local outfilename = filename:match('(.*)%.txt')..'.png'

local datatype = 'lines'
-- points for higher dimensions
-- if dim >= 3 then datatype = 'points' end

local cmds = {
	'set style data '..datatype,
	'set palette rgbformulae 33,13,10', --?
	'set output "'..outfilename..'"',
	'set format x "%e"',
	'set format y "%e"',
	'set format z "%e"'
}

local rows = {}
if useHistory then table.insert(rows, 1) end
for i=1,dim do
	table.insert(rows, i+1)
end

local i = 1
local found = false
for w in l:gmatch('[^\t]+') do
	if w == field then
		found = true
		table.insert(rows, i)
		break
	end
	i = i + 1
end
assert(found, "failed to find row")

local plotfunction = #rows <= 2 and 'plot' or 'splot'

-- if we're not in 3D already then add a color dimension
if #rows < 3 and useColor then table.insert(rows, rows[#rows]) end

local plotcmd = plotfunction..' "'..filename..'" using '..table.concat(rows, ':')
if useColor then plotcmd = plotcmd ..' with '..datatype..' palette' end
table.insert(cmds, plotcmd)

local gnuplot = 'droidplot'
local cmd = gnuplot.." -e '"..table.concat(cmds, '; ').."'"
print(cmd)
assert(os.execute(cmd) == 0)
