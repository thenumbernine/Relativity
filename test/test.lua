--[[
args
	setbase: 	set the current test results to the baseline.  notice that it does not run the tests and generate new results.
	run:		run tests
	col <col>	specify plot column
	plot:    	plot the results
	diff:		display warnings of what has changed
default if nothing is specified:
	run plot diff
--]]

local baselinePrefix = 'baseline_'
local plotColumn = 'K'
local diff = 'diff'
local plotLog = false		-- use log scaling
local plotHistory = true	-- whether to plot history or just the last state 


local args = {...}
local targetTestName = table.remove(args, 1)
print('targetTestName',targetTestName)
if #args == 0 then
	args = {'run', 'plot', 'diff'}
else
	local i = 1
	while i <= #args do
		if args[i] == 'col' then
			table.remove(args, i)
			plotColumn = table.remove(args, i)
		elseif args[i] == 'log' then
			table.remove(args, i)
			plotLog = true
		else
			i = i + 1
		end
	end
end

-- this one needs to match the makefile variable.
-- i might put that in a separate param file for both of these to read in
local installDir = '/data/local/bin/'	

local tests = {
	{
		name = 'gro_j0422_32_kerr-schild',
		dim = 1,
		iter = 100,
		res = 100,
		size = 4.1,
		args = 'kerr-schild 4.1',
	},
	{
		name = 'sagitarrius_a_star_kerr-schild',
		dim = 1,
		iter = 100,
		res = 100,
		size = 4.1e6,
		args = 'kerr-schild 4.1e6',
	},
	{
		name = 'binary_black_holes_brill-lindquist',
		dim = 1,
		iter = 100,
		res = 100,
		size = 4.1,
		args = 'brill-lindquist 2 -2 1 2 1', 
	},
	{
		name = 'gro_j0422_32_bowen-york',
		dim = 1,
		iter = 100,
		res = 100,
		size = 4.1,
		args = 'bowen-york 4.1 4.1',
	},
	{	-- spinning test
		name = 'sagitarrius_a_star_bowen-york',
		dim = 1,
		iter = 100,
		res = 100,
		size = 4.1e6,
		args = 'bowen-york 4.1e6 0 0 4.1e6',
	},
	{	-- boosted test
		name = 'sagitarrius_a_star_bowen-york_boosted',
		dim = 1,
		iter = 100,
		res = 100,
		size = 4.1e6,
		args = 'bowen-york 4.1e6 0 0 4.1e6 4.1e6 0 0',
	},
	{
		name = 'sagitarrius_a_star_bowen-york_2d',
		dim = 2,
		iter = 10,
		res = 20,
		size = 4.1,
		history = true,
		args = 'bowen-york 4.1e6 0 0 4.1e6',
	},
	{
		name = 'sagitarrius_a_star_bowen-york_2d_big',
		dim = 2,
		iter = 100,
		res = 100,
		size = 4.1,
		args = 'bowen-york 4.1e6 0 0 4.1e6',
	},

}

function io.fileexists(filename)
	local f = io.open(filename, 'r')
	if not f then return false end
	f:close()
	return true
end

local function exec(cmd)
	print(cmd)
	local result = os.execute(cmd)
	if result ~= 0 then return false, 'failed with error: '..result end
	return true
end

local function runCmd(targetTestName, cmd)
	local found = false
	for _,test in ipairs(tests) do
		if targetTestName == test.name
		or (targetTestName == 'all' and not test.excludeFromAll)
		then
			found = true
			local filename = test.name..'.txt'
			local basefile = baselinePrefix..filename
			if cmd == 'setbase' then
				assert(exec('cp '..filename..' '..basefile))
			end
			if cmd == 'run' then
				assert(exec(installDir..'relativity integrator rk4'
					..' filename '..filename
					..' dim '..test.dim
					..' iter '..test.iter
					..' res '..test.res
					..' size '..test.size
					..' allcols'
					..' history'
					..' '..test.args
				))
			end
			--[[
			plots are truly different than regression tests ...
			i only need a short number of iterations to run for regression tests before i see i've done something wrong
			but plots need more
			--]]
			if cmd == 'plot' then
				-- don't assert because gnuplot errors if the data is all zeroes
				exec('lua plot.lua '..filename
					..' '..(filename:match('(.*)%.txt'))..'.png'
					..' '..test.dim
					..' '..plotColumn
					..(plotLog and ' log' or '')
				)
				if plotHistory then
					exec('lua plot.lua '..filename
						..' '..(filename:match('(.*)%.txt'))..'_history.png'
						..' '..test.dim
						..' '..plotColumn
						..(plotLog and ' log' or '')
						..' history'
					)
				end
			end
			if cmd == 'diff' then
				if io.fileexists(basefile) then
					if not exec(diff..' '..basefile..' '..filename..' > /dev/null') then
						-- something went wrong? now we have to get our hands dirty...
						assert(exec('lua compare.lua '..filename..' '..basefile))
					end
				else
					io.stderr:write('baseline does not exist for '..filename..'\n')
					io.stderr:flush()
				end
			end
		end
	end
	if not found then
		io.stderr:write('failed to run test '..('%q'):format(targetTestName or '')..'\n')
		io.stderr:write('here are the valid test options:\n')
		for _,test in ipairs(tests) do
			io.stderr:write('\t'..('%q'):format(test.name)..'\n')
		end
		io.stderr:flush()
	end
end

for _,cmd in ipairs(args) do
	runCmd(targetTestName, cmd)
end

