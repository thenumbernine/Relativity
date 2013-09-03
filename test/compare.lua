--[[
baseline first, compared-to-file second

parses out table values and returns numerical info on the differing columns
--]]

local filenameA, filenameB = ...
assert(filenameA)
assert(filenameB)

local function stretch(range, value)
	if range[1] > value then range[1] = value end
	if range[2] < value then range[2] = value end
end

local ErrorColumn = {}
ErrorColumn.__index = ErrorColumn
setmetatable(ErrorColumn, {
	__call = function()
		local self = setmetatable({}, ErrorColumn)
		self.absARange = {math.huge, -math.huge}
		self.absBRange = {math.huge, -math.huge}
		return self
	end,
})
function ErrorColumn:process(na, nb)
	stretch(self.absARange, na)
	stretch(self.absBRange, nb)
	local relError = math.abs(nb - na) / math.abs(na)
end
function ErrorColumn:printStats(name)
	io.stderr:write(name..'\n')
	local absARange = self.absARange[2] - self.absARange[1]
	local absBRange = self.absBRange[2] - self.absBRange[1]
	local rangeRelErr = math.abs(absBRange - absARange) / math.abs(absARange)
	io.stderr:write('\tabs range in file A: '..self.absARange[1]..' to '..self.absARange[2]..'\n')
	io.stderr:write('\tabs range in file B: '..self.absBRange[1]..' to '..self.absBRange[2]..'\n')
	io.stderr:write('\trange relative error: '..rangeRelErr..'\n')
end

-- 1) make sure columns match
local fileA = assert(io.open(filenameA, 'r'))
local fileB = assert(io.open(filenameB, 'r'))

local function nextLines()
	return fileA:read('*l'), fileB:read('*l')
end

local function errorLines(lineA, lineB, msg)
	io.stderr:write(filenameA..': '..lineA..'\n')
	io.stderr:write(filenameB..': '..lineB..'\n')
	error(msg)
end

local headersA, headersB = nextLines()

if headersA ~= headersB then
	errorLines(headersA, headersB, "headers do not match!")
end

local function getColumns(line)
	local cols = {}
	for col in line:gmatch('%S+') do
		table.insert(cols, col)
	end
	return cols
end

local function getColumnNames(header)
	assert(header:sub(1,1) == '#', "expected a # at the beginning of the header")
	return getColumns(header:sub(2))
end

local headers = getColumnNames(headersA)	-- pick one, they're both the same

local errorColumns = {}
while true do
	if not fileA:read(0) or not fileB:read(0) then break end

	local lineA, lineB = nextLines()
	if lineA ~= lineB then
		-- split them up
		local colsA = getColumns(lineA)
		local colsB = getColumns(lineB)

		if #colsA ~= #headers then 
			errorLines(lineA, lineB, "column in file A has different number of elements than the header") 
		end
		if #colsB ~= #headers then 
			errorLines(lineA, lineB, "column in file B has different number of elements than the header") 
		end
	
		for i=1,#headers do
			if colsA[i] ~= colsB[i] then
				-- found a different element -- keep track of it
				-- keep track of ...
				-- 	range of relative error
				--  range of absolute error
				-- 	range of both values?
				--  coordinates of greatest relative error?
				local na = tonumber(colsA[i])
				local nb = tonumber(colsB[i])
				
				if not errorColumns[i] then errorColumns[i] = ErrorColumn() end
				errorColumns[i]:process(na, nb)
			end
		end
	end
end

if fileA:read(0) then
	io.stderr:write(filenameA..' has extra lines\n')
end
if fileB:read(0) then
	io.stderr:write(filenameB..' has extra lines\n')
end

for i,header in ipairs(headers) do
	local errorColumn = errorColumns[i]
	if errorColumn then
		errorColumn:printStats(headers[i])
	end
end
