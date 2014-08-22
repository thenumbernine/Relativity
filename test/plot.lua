#! /usr/bin/env lua
--[[
plotting script

plot.lua <filename.txt> <outfilename.png> <dim> <field> <flags>

flags:
	history
	log
	plotargs <plotargs>	-- arguments for the plot command
--]]

local args = {...}
local filename = table.remove(args, 1)
local outfilename = table.remove(args, 1)
local dim = table.remove(args, 1)
local field = table.remove(args, 1)
local plotargs = ''

dim = assert(tonumber(dim), 'expected a valid dimension')

local useHistory = false
local useLog = false
do
	local i = 1
	while i <= #args do
		local v = args[i]
		if v == 'history' then 
			useHistory = true 
		elseif v == 'log' then
			useLog = true
		elseif v == 'plotargs' then
			plotargs = table.remove(args, i+1)
		end 
		i = i + 1
	end
end
local useColor = true 

local f = assert(io.open(filename, 'rb'), "failed to open file")
local l = f:read('*l')
f:close()

assert(l:sub(1,1) == '#', "expected a comment on the first line")
l = l:sub(2)

local datatype = 'lines'
-- points for higher dimensions
-- if dim >= 3 then datatype = 'points' end

local cmds = {
	'set style data '..datatype,
	'set palette rgbformulae 33,13,10', --?
	'set terminal png size 800,600 enhanced font "Helvetica,12"',
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

if useLog then
	if #rows <= 1 then
		table.insert(cmds, 1, 'set log y')
	else
		table.insert(cmds, 1, 'set log z')
	end
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

local plotcmd = plotfunction..' '..plotargs..' "'..filename..'" using '..table.concat(rows, ':')
if useColor then plotcmd = plotcmd ..' with '..datatype..' palette' end
table.insert(cmds, plotcmd)

local gnuplot = 'gnuplot'
local cmd = gnuplot.." -e '"..table.concat(cmds, '; ').."'"
print(cmd)
assert(os.execute(cmd))
