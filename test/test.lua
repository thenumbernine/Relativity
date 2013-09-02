local baselinePrefix = 'baseline_'
local plotField = 'K'
local diff = 'diff'

-- this one needs to match the makefile variable.
-- i might put that in a separate param file for both of these to read in
local installDir = '/data/local/bin/'	


--[[
args
	setbase: 	set the current test results to the baseline.  notice that it does not run the tests and generate new results.
	run:		run tests
	plot:    	plot the results
	diff:		display warnings of what has changed
default if nothing is specified:
	run plot diff
--]]
local setBaseline = false
local testArgs = {}
for _,arg in ipairs{...} do 
	if arg == 'run' then testArgs.run = true end 
	if arg == 'setbase' then testArgs.setbase = true end 
	if arg == 'diff' then testArgs.warn = true end 
	if arg == 'plot' then testArgs.plot = true end 
end
-- if nothing was added then do the default:
local found = false
for k,v in pairs(testArgs) do found = true break end
if not found then
	testArgs.run = true
	testArgs.plot = true
	testArgs.diff = true
end

local tests = {
	{
		name = 'gro_j0422_32_kerr-schild',
		dim = 1,
		iter = 100,
		res = 100,
		history = true,
		size = 4.1,
		args = 'kerr-schild 4.1',
	},
	{
		name = 'sagitarrius_a_star_kerr-schild',
		dim = 1,
		iter = 100,
		res = 100,
		history = true,
		size = 4.1e6,
		args = 'kerr-schild 4.1e6',
	},
	{
		name = 'binary_black_holes_brill-lindquist',
		dim = 1,
		iter = 100,
		res = 100,
		history = true,
		size = 4.1,
		args = 'brill-lindquist 2 -2 1 2 1', 
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

function runTests(args)
	for _,test in ipairs(tests) do
		local filename = test.name..'.txt'
		local basefile = baselinePrefix..filename
		if args.setbase then
			assert(exec('cp '..filename..' '..basefile))
		end
		if args.run then
			assert(exec(installDir..'relativity integrator rk4'
				..' filename '..filename
				..' dim '..test.dim
				..' iter '..test.iter
				..' res '..test.res
				..' size '..test.size
				..(test.history and ' history' or '')
				..' '..test.args))
		end
		if args.plot then
			assert(exec('lua plot.lua '..filename
				..' '..test.dim
				..' '..plotField
				..(test.history and ' history' or '')))
		end
		if args.diff then
			if io.fileexists(basefile) then
				assert(exec(diff..' '..filename..' '..basefile))
			else
				io.stderr:write('baseline does not exist for '..filename..'\n')
				io.stderr:flush()
			end
		end
	end
end

runTests(testArgs)

