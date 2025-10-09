#!/bin/bash

# OpenFOAM Stabilized Windkessel Monitoring Script
LOG_FILE="logs/log.solver.final"

echo "======================================================"
echo "   STABILIZED WINDKESSEL MONITORING DASHBOARD"
echo "======================================================"
echo "Target endTime: 0.5s | Target timestep: 1e-5s"
echo ""

while true; do
    if [ ! -f "$LOG_FILE" ]; then
        echo "Waiting for log file to appear..."
        sleep 5
        continue
    fi
    
    # Clear screen for live updates
    clear
    echo "======================================================"
    echo "   STABILIZED WINDKESSEL MONITORING DASHBOARD"
    echo "======================================================"
    echo "$(date)"
    echo ""
    
    # Current simulation time and progress
    CURRENT_TIME=$(grep "Time =" $LOG_FILE | tail -1 | grep -o "Time = [0-9.e-]*" | cut -d' ' -f3)
    if [ ! -z "$CURRENT_TIME" ]; then
        PROGRESS=$(echo "scale=2; $CURRENT_TIME / 0.5 * 100" | bc 2>/dev/null || echo "0")
        echo "🎯 SIMULATION PROGRESS:"
        echo "   Current Time: ${CURRENT_TIME}s"
        echo "   Progress: ${PROGRESS}% of 0.5s target"
        echo ""
    fi
    
    # Timestep analysis
    echo "⏱️  TIMESTEP ANALYSIS:"
    RECENT_DT=$(grep "deltaT =" $LOG_FILE | tail -3)
    if [ ! -z "$RECENT_DT" ]; then
        echo "$RECENT_DT" | while read line; do
            DT_VAL=$(echo $line | cut -d'=' -f2 | tr -d ' ')
            echo "   $line"
        done
        LATEST_DT=$(echo "$RECENT_DT" | tail -1 | cut -d'=' -f2 | tr -d ' ')
        TARGET_RATIO=$(echo "scale=1; $LATEST_DT / 0.00001" | bc 2>/dev/null || echo "?")
        echo "   Ratio to target (1e-5): ${TARGET_RATIO}x"
    else
        echo "   No timestep data yet..."
    fi
    echo ""
    
    # Courant numbers
    echo "🌊 COURANT NUMBERS:"
    RECENT_CO=$(grep "Courant.*max" $LOG_FILE | tail -3)
    if [ ! -z "$RECENT_CO" ]; then
        echo "$RECENT_CO" | tail -1
        MAX_CO=$(echo "$RECENT_CO" | tail -1 | grep -o "max: [0-9.]*" | cut -d' ' -f2)
        if (( $(echo "$MAX_CO > 1.0" | bc -l 2>/dev/null || echo 0) )); then
            echo "   ⚠️  WARNING: Max Courant > 1.0"
        else
            echo "   ✅ Courant numbers stable"
        fi
    else
        echo "   No Courant data yet..."
    fi
    echo ""
    
    # Convergence info
    echo "🔄 CONVERGENCE STATUS:"
    PIMPLE_ITERS=$(grep "PIMPLE: Iteration" $LOG_FILE | tail -5 | wc -l)
    if [ $PIMPLE_ITERS -gt 0 ]; then
        LATEST_ITER=$(grep "PIMPLE: Iteration" $LOG_FILE | tail -1 | grep -o "Iteration [0-9]*" | cut -d' ' -f2)
        echo "   Latest PIMPLE iteration: $LATEST_ITER"
        if [ $LATEST_ITER -gt 20 ]; then
            echo "   ⚠️  High iteration count (may indicate slow convergence)"
        else
            echo "   ✅ Reasonable iteration count"
        fi
    fi
    
    # Continuity errors
    CONTINUITY=$(grep "time step continuity errors" $LOG_FILE | tail -1)
    if [ ! -z "$CONTINUITY" ]; then
        echo "   $CONTINUITY"
    fi
    echo ""
    
    # Performance estimates
    if [ ! -z "$CURRENT_TIME" ] && [ ! -z "$LATEST_DT" ]; then
        REMAINING_TIME=$(echo "scale=6; 0.5 - $CURRENT_TIME" | bc 2>/dev/null || echo "0.5")
        STEPS_REMAINING=$(echo "scale=0; $REMAINING_TIME / $LATEST_DT" | bc 2>/dev/null || echo "?")
        echo "📊 PERFORMANCE ESTIMATE:"
        echo "   Remaining physical time: ${REMAINING_TIME}s"
        echo "   Estimated steps remaining: $STEPS_REMAINING"
        echo ""
    fi
    
    # Stability indicators
    echo "🛡️  STABILIZATION STATUS:"
    echo "   Beta parameter: 1.0"
    echo "   Damping factor: 0.1"
    echo "   Backflow stabilization: ACTIVE"
    echo ""
    
    echo "Press Ctrl+C to stop monitoring..."
    echo "======================================================"
    
    sleep 10
done