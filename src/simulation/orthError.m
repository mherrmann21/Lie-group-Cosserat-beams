function err = orthError( R )
    %ORTHERROR Compute orthogonality error of a rotation matrix
    %   As done in [Lee08, p.75]
    %   R must be a 3x3 rotation matrix (element of SO3)

    arguments
        R (3,3)
    end

    err = norm( eye(3) - R' * R );
end

