function uiFibonaci(){
    let f = oneVirableFunctionParser(document.getElementById("fib-function").value)
    let epsilon = parseFloat(document.getElementById("fib-epsilon").value)
    let bottomInterval = parseFloat(document.getElementById("fib-bottomInterval").value)
    let topInterval = parseFloat(document.getElementById("fib-topInterval").value)
    
    let result = fibonaciSearch(f,bottomInterval,topInterval,epsilon)

    document.getElementById("fib-result").innerHTML = JSON.stringify(result,null,"\t")

}


function newtonConfig(){
    let data= {}
    data.epsilon = parseFloat(document.getElementById("newt-epsilon").value)
    data.startingPoint =parseFloat(document.getElementById("newt-start").value)
    data.numericalDerivatives = document.getElementById("newt-numericalDerivatives").checked
    return data
}

function uiNewton(){
    let config = newtonConfig()
    let func = oneVirableFunctionParser(document.getElementById("newt-function").value)

    let result

    if(config.numericalDerivatives){
        result = newtonMethod(func,config.startingPoint,config.epsilon)
    }else{
        let firstDeriv = oneVirableFunctionParser(document.getElementById("newt-firstDerivative").value)
        let secondDeriv = oneVirableFunctionParser(document.getElementById("newt-secondDerivative").value)
        result = generalNewtonAlgorythm(func,config.startingPoint,config.epsilon,firstDeriv,secondDeriv)
    }

    document.getElementById("newt-result").innerHTML = JSON.stringify(result,null,"\t")

}


function powellConfig(){
    let config = {}
    config.epsilon = parseFloat(document.getElementById("pow-epsilon").value)
    config.start = parseFloat(document.getElementById("pow-start").value)
    config.interval = parseFloat(document.getElementById("pow-interval").value)
    return config
}

function uiPowell(){
    let config = powellConfig()
    let f = oneVirableFunctionParser(document.getElementById("pow-function").value)

    let result = powellAlgorithm(f,config.start,config.interval,config.epsilon)

    document.getElementById("pow-result").innerHTML = JSON.stringify(result,null,"\t")
}

function getVector(inputId){
    return document.getElementById(inputId).value
        .replace(/\s/g,'')
        .split(",")
        .map(x=>parseFloat(x))
}

function uiGauss(){

    let f= functionParser(document.getElementById("gauss-function").value)
    let epsilon = parseFloat(document.getElementById("gauss-epsilon").value)

    let startingPoint = getVector("gauss-start")
    let result = gaussSeindel(f,startingPoint,epsilon)
    document.getElementById("gauss-result").innerHTML = JSON.stringify(result,null,"\t")

    if(startingPoint.length == 2){
        plotAlgorythm(f,result,1000)

    }

}

function uiFast(){
    let startingPoint = getVector("fast-start")
    let f = functionParser(document.getElementById("fast-function").value)
    let epsilon = parseFloat(document.getElementById("fast-epsilon").value)
    
    let result = fastestDesent(f,startingPoint,epsilon)
    
    document.getElementById("fast-result").innerHTML = JSON.stringify(result,null,"\t")
    console.log(result)

    if(startingPoint.length == 2){
        plotAlgorythm(f,result,1000)

    }


}

function uiPowc(){
    let startingPoint = getVector("powc-start")
    let f = functionParser(document.getElementById("powc-function").value)
    let epsilon = parseFloat(document.getElementById("powc-epsilon").value)
    
    let result = powellsCoupeldDirections(f,startingPoint,epsilon)
    
    document.getElementById("powc-result").innerHTML = JSON.stringify(result,null,"\t")
    console.log(result)

    if(startingPoint.length == 2){
        plotAlgorythm(f,result,1000)

    }

}


function uiCgrad(){
    let startingPoint = getVector("cgrad-start")
    let f = functionParser(document.getElementById("cgrad-function").value)
    let epsilon = parseFloat(document.getElementById("cgrad-epsilon").value)
    
    let result = coupeldGradient(f,startingPoint,epsilon)
    
    document.getElementById("cgrad-result").innerHTML = JSON.stringify(result,null,"\t")
    console.log(result)

    if(startingPoint.length == 2){
        plotAlgorythm(f,result,1000)

    }

}




function uiMin(){
    let config = newtonConfig()
    let func = functionParser(document.getElementById("min-function").value)


    let data= {}
    let epsilon = parseFloat(document.getElementById("min-epsilon").value)
    let startingPoint = getVector("min-start")
    let numericalDerivatives = document.getElementById("min-numericalDerivatives").checked

    let direction = getVector("min-direction")

    let result

    if(numericalDerivatives){
        result = minimumInDirection(func,direction,startingPoint,epsilon)
    }else{
        let firstDeriv = functionParser(document.getElementById("min-firstDerivative").value)
        let secondDeriv = functionParser(document.getElementById("min-secondDerivative").value)
        result = generalNewtonAlgorythm(func,direction,startingPoint,epsilon,firstDeriv,secondDeriv)
    }

    console.log(result)

    document.getElementById("min-result").innerHTML = JSON.stringify(result,null,"\t")



}

function uiGrad(){
     let func = functionParser(document.getElementById("grad-function").value)
     let point = getVector("grad-point")
     let delta = parseFloat(document.getElementById("grad-delta").value)

     document.getElementById("grad-result").innerHTML = JSON.stringify(gradient(func,point,delta),null,"\t")

}
