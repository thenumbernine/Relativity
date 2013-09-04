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
	local relMinErr = math.abs(self.absBRange[1] - self.absARange[1]) / math.abs(self.absARange[1])
	local relMaxErr = math.abs(self.absBRange[2] - self.absARange[2]) / math.abs(self.absARange[2])
	local relRangeErr = math.abs(absBRange - absARange) / math.abs(absARange)
	io.stderr:write('\tabs range in file A: '..self.absARange[1]..' to '..self.absARange[2]..'\n')
	io.stderr:write('\tabs range in file B: '..self.absBRange[1]..' to '..self.absBRange[2]..'\n')
	io.stderr:write('\trelative min error: '..relMinErr..'\n')
	io.stderr:write('\trelative max error: '..relMinErr..'\n')
	io.stderr:write('\trelative range error: '..relRangeErr..'\n')
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
	print("column names do not match!")
	io.stderr:write(filenameA..': '..headersA..'\n')
	io.stderr:write(filenameB..': '..headersB..'\n')
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

local colNamesA = getColumnNames(headersA)
local colNamesB = getColumnNames(headersB)

local colIndexForNameA = {}
local colIndexForNameB = {}
for i,name in ipairs(colNamesA) do colIndexForNameA[name] = i end
for i,name in ipairs(colNamesB) do colIndexForNameB[name] = i end

-- get common subset
local function getAllColNames(colNamesA, colNamesB)
	-- store as keys / set
	local allColNameKeys = {}
	for i,name in ipairs(colNamesA) do allColNameKeys[name] = allColNameKeys[name] and (allColNameKeys[name] + i) / 2 or i end
	for i,name in ipairs(colNamesB) do allColNameKeys[name] = allColNameKeys[name] and (allColNameKeys[name] + i) / 2 or i end
	-- convert to values / array
	-- maintain their index in the original so we can sort by that and somewhat preserve location (I'm sure there's a better way)
	local allColNames = {}
	for name,avgIndex in pairs(allColNameKeys) do
		table.insert(allColNames, {name=name, avgIndex=avgIndex})
	end
	-- sort by avg index
	table.sort(allColNames, function(a,b)
		return a.avgIndex < b.avgIndex
	end)
	-- extract names
	for i=1,#allColNames do allColNames[i] = allColNames[i].name end

	return allColNames
end
local allColNames = getAllColNames(colNamesA, colNamesB)
-- report any columns missing in any of the files
for _,name in ipairs(allColNames) do
	if not colIndexForNameA[name] then
		io.stderr:write("file A does not have column "..name.."\n")
	end
	if not colIndexForNameB[name] then
		io.stderr:write("file B does not have column "..name.."\n")
	end
end

local errorColumns = {}
while true do
	if not fileA:read(0) or not fileB:read(0) then break end

	local lineA, lineB = nextLines()
	if lineA ~= lineB then
		-- split them up
		local colsA = getColumns(lineA)
		local colsB = getColumns(lineB)

		for i,colName in ipairs(allColNames) do
			local indexA = colIndexForNameA[colName]
			local indexB = colIndexForNameB[colName]
			if indexA and indexB then
				if colsA[indexA] ~= colsB[indexB] then
					-- found a different element -- keep track of it
					-- keep track of ...
					-- 	range of relative error
					--  range of absolute error
					-- 	range of both values?
					--  coordinates of greatest relative error?
					local na = tonumber(colsA[indexA])
					local nb = tonumber(colsB[indexB])
					
					if not errorColumns[i] then errorColumns[i] = ErrorColumn() end
					errorColumns[i]:process(na, nb)
				end
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

for i,colName in ipairs(allColNames) do
	local errorColumn = errorColumns[i]
	if errorColumn then
		errorColumn:printStats(allColNames[i])
	end
end
