pipeline {
    agent any

    environment {
        NEW_VERSION = '0.0.1'
    }
    
    stages {
        stage('build') {
            steps{
                echo 'Build'
                echo "Building version ${NEW_VERSION}"
            }
        }
        stage('test') {
            steps {
                echo 'Test'
            }
        }
        stage('deploy') {
            steps {
                echo 'Deploy'
            }
        }
    }

    post {
        always {

        }
        failure {

        }
    }
}