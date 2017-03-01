
function fibonaciNumber(number){
        if(number == 0) return 1
        if(number == 1) return 1

        let a = 1
        let b = 1
        let c


        for(let i = 0; i < number-1;i++){
            c = a+b
            b=a
            a=c
        }

        return c
}
function fibonaciSearch(functionToMinimize,bottomSearchInterval,topSearchInterval,epsilon){
    function findFibonaciIndex(x){
        let i = 0 ;
        while (fibonaciNumber(i)<x) i++
        return i
    }


    let data = {
        iterations:[],
    }
    let k = findFibonaciIndex((topSearchInterval - bottomSearchInterval)/epsilon)


    let x1
    let x2
    for (let i = 2; i <= k; i++){
        
        let delta = Math.abs(topSearchInterval - bottomSearchInterval)
        let alpha = fibonaciNumber(k-i)/fibonaciNumber(k-i+2)
        
         x1 = bottomSearchInterval+alpha*delta
         x2 = bottomSearchInterval + (1-alpha)*delta

         data.iterations.push({
            x1:x1,
            x2:x2,
            delta:delta,
            alpha:alpha,
            fx1:functionToMinimize(x1),
            fx2:functionToMinimize(x2),
            a:bottomSearchInterval,
            b:topSearchInterval,
            fib:[
                fibonaciNumber(k-i),fibonaciNumber(k-i+2)
            ]
         })




        if(functionToMinimize(x1)<functionToMinimize(x2)){
            topSearchInterval = x2
        }else {
            bottomSearchInterval = x1;
        }
    }
    data.result = (x1+x2)/2

    return data
}


function goldenRatioSearch(functionToMinimize,bottomSearchInterval,topSearchInterval){

    const goldenRation = 1.61803398875
    const alpha = 1-(1/goldenRation)

    let k = 6
    let x1
    let x2
    for(let i = 0;i<k;i++){
        let delta = Math.abs(topSearchInterval - bottomSearchInterval)
            x2 = topSearchInterval-delta*alpha
            x1 = bottomSearchInterval+delta*alpha

        if(functionToMinimize([x1])<=functionToMinimize([x2])){
                topSearchInterval = x2
            }else {
                bottomSearchInterval = x1;
            }
        }
        return (x1+x2)/2
}

function powellAlgorithm(functionToMinimize,startingPoint,interval,epsilon){
    let data = {
        iterations:[]
    }
    return _powellAlgorithm(functionToMinimize,startingPoint,interval,epsilon,data)

    function goRigth(startingPoint,interval,functionToMinimize,iteration){

        let leftPoint = startingPoint
        let rigthPoint = startingPoint + interval
        
        let leftPointValue = functionToMinimize(leftPoint)
        let rigthPointValue = functionToMinimize(rigthPoint)

        iteration.pointsSearch = [];
        iteration.pointsSearch.push({
            leftPoint:leftPoint,
            rigthPoint:rigthPoint,
            t:interval
        })

        while(leftPointValue > rigthPointValue){
            leftPoint = rigthPoint
            interval*=2
            rigthPoint +=interval 
            leftPointValue = functionToMinimize(leftPoint)
            rigthPointValue = functionToMinimize(rigthPoint)
            
            iteration.pointsSearch.push({
                leftPoint:leftPoint,
                rigthPoint:rigthPoint,
                t:interval
            })

        }
        rigthPoint = leftPoint
        interval/=2
        leftPoint -= interval
        let middlPoint = leftPoint + 0.5 * interval

        iteration.pointsSearch.push({
            leftPoint:leftPoint,
            rigthPoint:rigthPoint,
            t:interval
        })


        return {
            leftPoint: leftPoint,
            middlPoint : middlPoint,
            rigthPoint : rigthPoint
        }
        
    }
    function goLeft(startingPoint,interval,functionToMinimize,iteration){
        let rigthPoint = startingPoint
        let leftPoint = startingPoint - interval

        let leftPointValue = functionToMinimize(leftPoint)
        let rigthPointValue = functionToMinimize(rigthPoint)
                
        iteration.pointsSearch = [];
        iteration.pointsSearch.push({
                leftPoint:leftPoint,
                rigthPoint:rigthPoint,
                t:interval
            })

        while(rigthPointValue>leftPointValue){
            rigthPoint = leftPoint
            interval*=2
            leftPoint -= interval
            leftPointValue = functionToMinimize(leftPoint)
            rigthPointValue = functionToMinimize(rigthPoint)

            iteration.pointsSearch.push({
                leftPoint:leftPoint,
                rigthPoint:rigthPoint,
                t:interval
            })
        }        

        leftPoint = rigthPoint
        interval/=2
        rigthPoint +=interval
        let middlPoint = rigthPoint - 0.5 * interval

        iteration.pointsSearch.push({
                leftPoint:leftPoint,
                rigthPoint:rigthPoint,
                t:interval
            })


        return {
            leftPoint: leftPoint,
            middlPoint : middlPoint,
            rigthPoint : rigthPoint
        }

    }
    function _powellAlgorithm(functionToMinimize,startingPoint,interval,epsilon,data){

        let baseInterval = interval
        let x1 = startingPoint
        let x3 = x1+interval

        let y1 = functionToMinimize(x1)
        let y3 = functionToMinimize(x3)
        let result

        let iteration = {

        }

        if(y1>y3){
            result = goRigth(startingPoint,interval,functionToMinimize,iteration)
        }else{
            result = goLeft(startingPoint,interval,functionToMinimize,iteration)
        }

        let mat = [
                        [1, result.leftPoint, result.leftPoint*result.leftPoint],
                        [1, result.middlPoint, result.middlPoint*result.middlPoint],
                        [1, result.rigthPoint, result.rigthPoint*result.rigthPoint]
                    ]
        let matX = [
                        functionToMinimize(result.leftPoint),
                        functionToMinimize(result.middlPoint),
                        functionToMinimize(result.rigthPoint)
                    ]

        let interpolationParabola  = gauss(mat,matX)
        let parabolaA = interpolationParabola[2]
        let parabolaB = interpolationParabola[1]
        let parabolaC = interpolationParabola[0]
        
        let parabolaMinimumX = -parabolaB/(2*parabolaA)
        let parabolaMinumumY = parabolaA*parabolaMinimumX*parabolaMinimumX+ parabolaB*parabolaMinimumX + parabolaC

        let q = (parabolaMinumumY-functionToMinimize(parabolaMinimumX))/functionToMinimize(parabolaMinimumX)

        iteration.parabolaA = parabolaA
        iteration.parabolaB = parabolaB
        iteration.parabolaC = parabolaC
        iteration["parabola minimum x"] = parabolaMinimumX

        data.iterations.push(iteration);


        if(Math.abs(q) > epsilon){
            return _powellAlgorithm(functionToMinimize,parabolaMinimumX,interval/2,epsilon,data)
        }else
            data.result = parabolaMinimumX
            return data
    }


}


function generalNewtonAlgorythm(functionToMinimize,startingPoint,epsilon,derivative1,derivative2){
    
    let data = {
        iterations:[]

    }

    let currentPoint = startingPoint

    let previusPoint

    do{
        previusPoint = currentPoint


        currentPoint = currentPoint - (derivative1(currentPoint,epsilon)/
                                        derivative2(currentPoint,epsilon))
        data.iterations.push({
            x0:currentPoint,
            poch1 : derivative1(currentPoint,epsilon),
            poch2 : derivative2(currentPoint,epsilon)
        })

    }while(Math.abs(currentPoint-previusPoint)>epsilon)
    data.result = currentPoint

    return data

}


function newtonMethod(functionToMinimize,startingPoint,epsilon){
    let firstDerit = (x,delta)=>oneVirableFirstDerivative(functionToMinimize,x,delta)
    let secondDerit= (x,delta)=>oneVirableSecondDerivative(functionToMinimize,x,delta)

    return generalNewtonAlgorythm(functionToMinimize,startingPoint,epsilon,firstDerit,secondDerit)

}



function gaussSeindel(functionToMinimize,startingPoint,epsilon){

    let data = {
        epsilon:epsilon,
        points:[],
        iterations:[]
    }
    let numberOfDimension = startingPoint.length

    let singleVirableMinimization = _.curry(newtonMethod)(_,_,epsilon)

    let currentPos = startingPoint
    data.points.push(currentPos.slice())
    let previousCurrentPos
    let directions = identityMatrix(numberOfDimension)

    do{
        let iterationData = {
            dir:[]
        }

        previousCurrentPos = currentPos.slice()

        directions.forEach(direction=>{
            currentPos = minimumInDirection(functionToMinimize,direction,currentPos,epsilon).result
            data.points.push(currentPos.slice())
            iterationData.dir.push({
                position:currentPos.slice(),
                inDirection:direction.slice()
            })
        })
        data.iterations.push(iterationData)

    }while (distanceBetweenPoints(currentPos,previousCurrentPos)>epsilon)
    data.result = currentPos
    return data

}

function virableIndexArrayParameter(f,array,index){
    return function(x){
        array[index] = x
        return f(array)
    }
}



function fastestDesent (functionToMinimize,startingPoint,epsilon){

    let data = {
        points:[startingPoint],
        iterations:[]
    }

    let currentPos = startingPoint
    let previousCurrentPos
    do{
        let grad = gradient(functionToMinimize,currentPos,epsilon/2)      
        previousCurrentPos = currentPos.slice()

        currentPos = minimumInDirection(functionToMinimize,grad,currentPos,epsilon).result
        data.points.push(currentPos)
        data.iterations.push({
            position:currentPos.slice(),
            inDirection:grad
        })


    }while(distanceBetweenPoints(currentPos,previousCurrentPos)>epsilon)
    data.result = currentPos
    return data

}

function directionalSubfunction(func,direction,startingPoint){
    return (x)=>{
        return func(direction.map((dir,index)=>startingPoint[index]+dir*x))
    }
}

function minimumInDirection(func,dir,pos,epsilon){

    dir = normalizeVector(dir)

    const directionalFunction = directionalSubfunction(func,dir,pos)
    
    const firstDerit = (x,delta)=>oneVirableFirstDerivative(directionalFunction,x,delta)
    const secondDerit= (x,delta)=>oneVirableSecondDerivative(directionalFunction,x,delta)

    return generalMinimumInDirection(func,dir,pos,epsilon,firstDerit,secondDerit)

}

function generalMinimumInDirection(func,dir,pos,epsilon,firstDerit,secondDerit){
    
    let singleVirableMinimization = _.curry(generalNewtonAlgorythm)(_,_,epsilon,firstDerit,secondDerit)
    let directionalFunction = directionalSubfunction(func,dir,pos)
    let result = singleVirableMinimization(directionalFunction,0)

    return {
        minimization: result,
        result:pos.map((x,index)=>x+result.result*dir[index])

    }


}

function powellsCoupeldDirections(functionToMinimize,startingPoint,epsilon){

    let data = {
        points:[],
        iterations:[]
    }
    data.points.push(startingPoint.slice())
    
    let currentPos = startingPoint
    let singleVirableMinimization = _.curry(newtonMethod)(_,_,epsilon)
    let previousCurrentPos
    let numberOfVirables = startingPoint.length
    
    let directions = identityMatrix(numberOfVirables)

    do{
        let iteration={
            it:[]
        }

       let iterationPoints = []

       previousCurrentPos = currentPos.slice()
       currentPos = minimumInDirection(functionToMinimize,directions[numberOfVirables-1],currentPos,epsilon).result
       iterationPoints.push(currentPos.slice())
       data.points.push(currentPos.slice())

    for(let i = 0 ; i < numberOfVirables;i++){
        currentPos = minimumInDirection(functionToMinimize,directions[i],currentPos,epsilon).result
        data.points.push(currentPos.slice())
        iterationPoints.push(currentPos.slice())
        iteration.it.push({
            position:currentPos.slice(),
            inDirection:directions[i].slice()
        })
    }
    

    for(let i = 0;i<numberOfVirables-1;i++){
        directions[i] = directions[i+1]

       data.points.push(currentPos.slice())
       }
       
       directions[numberOfVirables-1] = normalizeVector(substractVectors(iterationPoints[numberOfVirables],iterationPoints[0]))

       iteration.directions = directions.map(x=>x.slice())
       data.iterations.push(iteration)

    }while(distanceBetweenPoints(currentPos,previousCurrentPos)>epsilon)

    return data

}

function coupeldGradient(functionToMinimize,startingPoint,epsilon){

    function nextDirection(direction,previousDirection){
        let beta = vectorLengthSquared(direction)/vectorLengthSquared(previousDirection)
        let dir = direction.map((x,index)=>-x+beta*previousDirection[index])
        return {
                    direction : dir,
                    beta:beta
                }
    }

    let data = {
        points:[],
        iterations:[]
    }
    data.points.push(startingPoint)

    let localGradient = (x)=>gradient(functionToMinimize,x,epsilon)

    let currentPos = startingPoint
    let previousCurrentPos

    let currentDirection = localGradient(currentPos).map(x=>-x)

    let next
    let beta = NaN
    do{
       
        
        previousCurrentPos = currentPos.slice() 
        currentPos = minimumInDirection(functionToMinimize,currentDirection,currentPos,epsilon).result
        data.points.push(currentPos)

        data.iterations.push({
            position:currentPos.slice(),
            inDirection:currentDirection.slice(),
            beta:beta
        })

        next = nextDirection(localGradient(currentPos),currentDirection,epsilon)
        currentDirection = next.direction
        beta = next.beta

    }while(distanceBetweenPoints(currentPos,previousCurrentPos)>epsilon)

    return data

}


function identityMatrix(n){
        return Array.apply(null, new Array(n))
        .map(function (x, i, xs) {
            return xs.map(function (_, k) {
                return i === k ? 1 : 0;
            })
        });
}

function distanceBetweenPoints(point1,point2){
    return Math.sqrt(point1.reduce((pv,cv,index)=>pv+ Math.pow(cv-point2[index],2),0))
}


function firstDerivative(func,point,virableIndex,delta){
    let x = point[virableIndex]
    let f = virableIndexArrayParameter(func,point,virableIndex)
    return oneVirableFirstDerivative(f,x,delta)
}

function secondDerivative(func,point,virableIndex,delta){
    let x = point[virableIndex]
    let f = virableIndexArrayParameter(func,point,virableIndex)
    return oneVirableSecondDerivative(f,x,delta)
}

function oneVirableFirstDerivative(f,x,delta){
    return (-(f(x + 2 * delta)) + 4 * (f(x + delta)) - 3 * f(x)) / (2 * delta)
}
function oneVirableSecondDerivative(f,x,delta){
    return (-f(x + 3 * delta) + 4 * f(x + 2 * delta) - 5 * f(x + delta) + 2 * f(x)) /(delta*delta)
}



function gradient(func,point,delta){
    return point.map((value,index)=>firstDerivative(func,point,index,delta))
}

function vectorLength(vector){
      return  Math.sqrt(vectorLengthSquared(vector))
}
function vectorLengthSquared(vector){
    return vector.reduce((pv, cv) => pv+(cv*cv), 0)
}

function normalizeVector(vector){
    let length = vectorLength(vector)
    return vector.map(x=>x/length)
}


function oneVirableFunctionParser(input){
    let f = functionParser(input)
    return (x)=>f([x])

}

function functionParser(input){
    input = input.replace(/\s/g,'')

    return createNode(input);

    function createNode(input){
        return node(stripOutsideParenthesis(input))

    }

    function stripOutsideParenthesis(input){
        if(input[0] == "(" && input[input.length-1] == ")" &&
            corespondingParenthesisIndex(input,0) === input.length-1){
            return stripOutsideParenthesis(input.substring(1,input.length-1))
        }else return input
    }


    

    function node(input){
    const predefineFunctions = {
        PI: () => Math.PI,
        E : ()=> Math.E,
        sin: x => Math.sin(x),
        cos: x => Math.cos(x),
        log: x => Math.log10(x),
        tan: x => Math.tan(x),
        abs: x => Math.abs(x),
        exp : x=> Math.exp(x),
        ln:  x => Math.log(x),

     }


        const functionRegex = /^(\w+)\((.*)\)$/
        const virableRegex = /^x(\d)$/

        if(!isNaN(input)){
            return (x) => parseFloat(input) //node is a leaf
        }else if (virableRegex.test(input)){
            let virableIndex = virableRegex.exec(input)[1]
            return (x)=>x[virableIndex];
        }else if (isFunction(input)){
            const groups = functionRegex.exec(input)
            const functionName = groups[1]
            const functionBody = groups[2]
            if(!predefineFunctions[functionName]) throw "function " + functionName + " does not exists"
            return (x)=> predefineFunctions[functionName]( createNode(functionBody)(x))
        }
        else{
            const result = splitToTwoExpresions(input)
            const children1 = createNode(result.leftExpresion)
            const children2 = createNode(result.rigthExpresion)
            const operation = result.operation
            return (x)=> operation(children1(x),children2(x))
        }
        
    }

    function isFunction(input){
        

        if(!/\w/.test(input[0])) return false
        for(let i = 1 ; i < input.length;i++){
            if(input[i]==="("){
                if(corespondingParenthesisIndex(input,i)==input.length-1) return true;
                else return false
            }
            if(!/\w/.test(input[i])) return false
        }
        return false;
    }

    function corespondingParenthesisIndex(input,bracketIndex){
        let bracketCount = 0;
        for(let i = bracketIndex;i<input.length;i++){
            if(input[i]==="(") bracketCount++;
            else if(input[i]===")") bracketCount--;
            if(bracketCount == 0) return i
        }
        throw "No coresponding bracket "
    }

    function splitToTwoExpresions(input){

        function getOperation(char){
            const operations = {
                "+": (x,y) => x+y,
                "-": (x,y) => x-y,
                "*": (x,y) => x*y,
                "/": (x,y) => x/y,
                "^": (x,y) => Math.pow(x,y),                                               
            }
            return operations[char]
        }

        
        function findNextOperator(input,startIndex){
            for(let i = startIndex;i<input.length;i++){
                if(input[i]==="("){ 
                    i = corespondingParenthesisIndex(input,i)+1
                }
                if(isOperator(input[i])) return i
            }
            return -1;
        }

        let indexToSplit = findNextOperator(input,0)
        let nextOperatorIndex = findNextOperator(input,indexToSplit+1)
        while( nextOperatorIndex !==-1){
            if(compareOperaionOrder(input[indexToSplit],input[nextOperatorIndex])>=0){
                indexToSplit = nextOperatorIndex
            }
            nextOperatorIndex = findNextOperator(input,nextOperatorIndex+1)

        }


        if (indexToSplit ===-1) throw "error while parsing"

        let leftExpresion = input.substring(0,indexToSplit)
        let rigthExpresion = input.substring(indexToSplit+1)

        return {
                    operation : getOperation(input[indexToSplit]),
                    leftExpresion: leftExpresion,
                    rigthExpresion :rigthExpresion
                }
    }


    function isOperator(char){
        return isLowerOrderOperation(char) || isMiddleOrderOperation(char) || isHigherierOrderOperation(char)
    }

    function isLowerOrderOperation(char){
        if(char ==='+' || char == '-')
        return char ==='+' || char == '-'
    }
    function isMiddleOrderOperation(char){
        return char ==='*' || char === '/'
    }

    function isHigherierOrderOperation(char){
        return char ==='^'
    }

    function compareOperaionOrder(a,b){
        const orderTable = {
            "+":0,
            "-":0,
            "*":1,
            "/":1,
            "^":2,
        }
        return orderTable[a]-orderTable[b];
    }


}

function substractVectors(v1,v2){
    return v1.map((x,index)=>x-v2[index])
}


function plotAlgorythm(func,data,numberOfPoints){
    function minValue(arr,index){
        let min = arr[0][index]
        for(let i = 1; i < arr.length;i++){
            if(arr[i][index]<min) 
               min =  arr[i][index]
        }
        return min
    }
    function maxValue(arr,index){
        let max = arr[0][index]
        for(let i = 1; i < arr.length;i++){
            if(arr[i][index]>max) 
               max = arr[i][index]
        }
        return max
    }
    function funcSurface(func,startX,startY,length,epsilon){
        let surface = {
            x:[],
            y:[],
            z:[],
            i:[],
            j:[],
            k:[],
            type:"mesh3d",
            lighting:{
               // fresnel:0.9,
               // facenormalsepsilon:1e-3,
               // roughness:0,
               // specular:0,
                ambient:0.2,
                diffuse:1
            },
            opacity:0.9


        }
        for(let i = startX;i< startX+length;i+=epsilon){
           let row = []
            for(let j = startY;j< startY+length;j+=epsilon){
                surface.x.push(j)
                surface.y.push(i)
                surface.z.push(func([j,i]))

            }

        }
        let rowLength = Math.ceil(length/epsilon)

        for(let i = 0; i < rowLength - 1;i++){
            for(let j = 0; j < rowLength - 1;j++){

                surface.j.push(i*(rowLength)+j)
                surface.k.push(i*(rowLength)+j+1)
                surface.i.push((i+1)*rowLength+j) 

                surface.i.push(i*(rowLength)+j+1)
                surface.k.push((i+1)*rowLength+j) 
                surface.j.push((i+1)*rowLength+j+1) 
            }
        }



        return surface
    }

    function scaterData(func,points){
        let scaterGraph = {
            x:[],
            y:[],
            z:[],
            type:"scatter3d",
            
        }

        for(let i = 0 ;i< points.length;i++){
            scaterGraph.x.push(points[i][0])
            scaterGraph.y.push(points[i][1])
            scaterGraph.z.push(func(points[i]))
        }
        return scaterGraph
    }


    let minX = minValue(data.points,0)
    let minY = minValue(data.points,1)
    let maxX = maxValue(data.points,0)
    let maxY = maxValue(data.points,1)

    let lengthX = Math.abs(minX - maxX)
    let lengthY = Math.abs(minY - maxY)
    
    let length = Math.max(lengthX,lengthY)

    let margin = length/2
    let totalLength = length+(margin*2)

    let epsilon = totalLength/Math.sqrt(numberOfPoints)

    let surface = funcSurface(func,
                    minY-margin,
                    minX-margin,
                    length+margin*2,    
                    epsilon)

    let dat = []
    dat.push(surface)
    dat.push(scaterData(func,data.points))

    let layout = {
       width: 900,
       height: 900,
    }

    Plotly.purge("plotly")
    Plotly.newPlot('plotly',dat,layout)

}








































/**
 * Gaussian elimination
 * @param  array A matrix
 * @param  array x vector
 * @return array x solution vector
 */
function gauss(A, x) {
    var abs = Math.abs;

    function array_fill(i, n, v) {
        var a = [];
        for (; i < n; i++) {
            a.push(v);
        }
        return a;
    }

    var i, k, j;

    // Just make a single matrix
    for (i=0; i < A.length; i++) { 
        A[i].push(x[i]);
    }
    var n = A.length;

    for (i=0; i < n; i++) { 
        // Search for maximum in this column
        var maxEl = abs(A[i][i]),
            maxRow = i;
        for (k=i+1; k < n; k++) { 
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }


        // Swap maximum row with current row (column by column)
        for (k=i; k < n+1; k++) { 
            var tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (k=i+1; k < n; k++) { 
            var c = -A[k][i]/A[i][i];
            for (j=i; j < n+1; j++) { 
                if (i===j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    x = array_fill(0, n, 0);
    for (i=n-1; i > -1; i--) { 
        x[i] = A[i][n]/A[i][i];
        for (k=i-1; k > -1; k--) { 
            A[k][n] -= A[k][i] * x[i];
        }
    }

    return x;
}
